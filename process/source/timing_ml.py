#!/usr/bin/env python3
"""Differentiable TIMING scintillator calibration/estimator study.

The model learns per-channel offsets and event-by-event channel weights for
TIMING0 and TIMING1 by minimizing the robust width of

    timing0_estimator - timing1_estimator

No external target time is used. The gauge is fixed by keeping offset0[0] = 0.
"""

import argparse
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


EO2DO = np.array(
    [
        22, 20, 18, 16, 24, 26, 28, 30,
        25, 27, 29, 31, 23, 21, 19, 17,
        9, 11, 13, 15, 7, 5, 3, 1,
        6, 4, 2, 0, 8, 10, 12, 14,
    ],
    dtype=np.int64,
)


def do_xy():
    x = EO2DO % 4
    y = EO2DO // 4
    return x.astype(np.float32), y.astype(np.float32)


def load_npz(path):
    data = np.load(path)
    return (
        data["t0"].astype(np.float64),
        data["m0"].astype(bool),
        data["t1"].astype(np.float64),
        data["m1"].astype(bool),
    )


def load_root_with_uproot(path):
    try:
        import uproot
        import awkward as ak
    except ImportError as exc:
        raise RuntimeError(
            "ROOT input requires uproot/awkward. Install uproot or pass --input-npz."
        ) from exc

    tree = uproot.open(path)["frames"]
    arrays = tree.arrays(
        [
            "ntiminghits",
            "timing_frame_start",
            "timing_frame_nhits",
            "timing_fifo",
            "timing_column",
            "timing_pixel",
            "timing_time",
        ],
        library="ak",
    )

    t0 = []
    t1 = []
    m0 = []
    m1 = []

    for entry in range(len(arrays["ntiminghits"])):
        nframes = len(arrays["timing_frame_start"][entry])
        starts = arrays["timing_frame_start"][entry]
        nhits = arrays["timing_frame_nhits"][entry]
        fifo = arrays["timing_fifo"][entry]
        column = arrays["timing_column"][entry]
        pixel = arrays["timing_pixel"][entry]
        time = arrays["timing_time"][entry]

        for iframe in range(nframes):
            row_t = np.zeros((2, 32), dtype=np.float64)
            row_m = np.zeros((2, 32), dtype=bool)
            first = int(starts[iframe])
            last = first + int(nhits[iframe])
            for i in range(first, last):
                f = int(fifo[i])
                det = 0 if 0 <= f <= 3 else 1 if 4 <= f <= 7 else -1
                if det < 0:
                    continue
                ch = int(pixel[i]) + 4 * int(column[i])
                if ch < 0 or ch >= 32:
                    continue
                value = float(time[i])
                if not row_m[det, ch] or value < row_t[det, ch]:
                    row_t[det, ch] = value
                    row_m[det, ch] = True
            t0.append(row_t[0])
            m0.append(row_m[0])
            t1.append(row_t[1])
            m1.append(row_m[1])

    return (
        np.asarray(t0, dtype=np.float64),
        np.asarray(m0, dtype=bool),
        np.asarray(t1, dtype=np.float64),
        np.asarray(m1, dtype=bool),
    )


def normalize_event_times(t0, m0, t1, m1):
    """Remove the large common event time before float32 training.

    ROOT hit times can be O(1e8) clocks, while the useful timing structure is
    O(1e-1) clocks.  Keeping the absolute timestamp in float32 destroys the
    fine timing information, so each event is shifted by its earliest TIMING
    hit before building tensors.  This does not change any delta time.
    """
    big0 = np.where(m0, t0, np.inf)
    big1 = np.where(m1, t1, np.inf)
    base = np.minimum(big0.min(axis=1), big1.min(axis=1))
    t0 = np.where(m0, t0 - base[:, None], 0.0)
    t1 = np.where(m1, t1 - base[:, None], 0.0)
    return t0.astype(np.float32), t1.astype(np.float32)


class TimingEstimator(nn.Module):
    def __init__(self, hidden=32, temperature=0.12, initial_offset0=None, initial_offset1=None):
        super().__init__()
        if initial_offset0 is None:
            initial_offset0 = torch.zeros(32)
        if initial_offset1 is None:
            initial_offset1 = torch.zeros(32)
        self.register_buffer("initial_offset0", initial_offset0.float())
        self.register_buffer("initial_offset1", initial_offset1.float())
        self.delta_offset0_free = nn.Parameter(torch.zeros(31))
        self.delta_offset1 = nn.Parameter(torch.zeros(32))
        self.temperature = temperature

        x, y = do_xy()
        geom = np.stack([x / 3.0, y / 7.0], axis=1)
        self.register_buffer("geom", torch.tensor(geom, dtype=torch.float32))

        self.score = nn.Sequential(
            nn.Linear(6, hidden),
            nn.Tanh(),
            nn.Linear(hidden, hidden),
            nn.Tanh(),
            nn.Linear(hidden, 1),
        )

    def offsets(self):
        delta0, delta1 = self.offset_deltas()
        return self.initial_offset0 + delta0, self.initial_offset1 + delta1

    def offset_deltas(self):
        delta0 = torch.cat(
            [
                torch.zeros(1, device=self.delta_offset0_free.device, dtype=self.delta_offset0_free.dtype),
                self.delta_offset0_free,
            ]
        )
        return delta0, self.delta_offset1

    def estimate_one(self, t, mask, offset):
        tc = t - offset
        big = torch.full_like(tc, 1.0e9)
        first = torch.min(torch.where(mask, tc, big), dim=1).values
        dt = tc - first[:, None]

        first_index = torch.argmin(torch.where(mask, tc, big), dim=1)
        first_geom = self.geom[first_index]
        geom = self.geom[None, :, :].expand(t.shape[0], -1, -1)
        rel_geom = geom - first_geom[:, None, :]

        features = torch.cat(
            [
                dt[:, :, None],
                torch.abs(dt[:, :, None]),
                geom,
                rel_geom,
            ],
            dim=2,
        )
        logits = self.score(features).squeeze(-1) / self.temperature
        logits = torch.where(mask, logits, torch.full_like(logits, -1.0e9))
        weights = torch.softmax(logits, dim=1)
        estimate = torch.sum(weights * tc, dim=1)
        neff = 1.0 / torch.sum(weights * weights, dim=1)
        return estimate, weights, neff

    def forward(self, t0, m0, t1, m1):
        offset0, offset1 = self.offsets()
        e0, w0, neff0 = self.estimate_one(t0, m0, offset0)
        e1, w1, neff1 = self.estimate_one(t1, m1, offset1)
        return e0, e1, w0, w1, neff0, neff1


def robust_center(x):
    return torch.quantile(x.detach(), 0.5)


def robust_loss(delta, huber=0.08):
    center = robust_center(delta)
    return torch.nn.functional.huber_loss(delta - center, torch.zeros_like(delta), delta=huber)


def evaluate(model, loader, device):
    model.eval()
    deltas = []
    neff0 = []
    neff1 = []
    with torch.no_grad():
        for t0, m0, t1, m1 in loader:
            t0 = t0.to(device)
            t1 = t1.to(device)
            m0 = m0.to(device)
            m1 = m1.to(device)
            e0, e1, _, _, n0, n1 = model(t0, m0, t1, m1)
            deltas.append((e0 - e1).cpu())
            neff0.append(n0.cpu())
            neff1.append(n1.cpu())
    delta = torch.cat(deltas).numpy()
    n0 = torch.cat(neff0).numpy()
    n1 = torch.cat(neff1).numpy()
    return {
        "entries": len(delta),
        "mean": float(delta.mean()),
        "rms": float(delta.std()),
        "median": float(np.median(delta)),
        "q16": float(np.quantile(delta, 0.16)),
        "q84": float(np.quantile(delta, 0.84)),
        "neff0": float(n0.mean()),
        "neff1": float(n1.mean()),
    }


def collect_diagnostics(model, loader, device, nbins=160, dt_range=(-2.0, 2.0)):
    model.eval()

    weight_sum = np.zeros((2, 32), dtype=np.float64)
    weight_by_first_sum = np.zeros((2, 32, 32), dtype=np.float64)
    first_count = np.zeros((2, 32), dtype=np.int64)
    neff_values = [[], []]
    dt_edges = np.linspace(dt_range[0], dt_range[1], nbins + 1, dtype=np.float64)
    dt_weight_sum = np.zeros((2, nbins), dtype=np.float64)
    dt_weight_count = np.zeros((2, nbins), dtype=np.int64)
    delta_timing = []

    with torch.no_grad():
        for t0, m0, t1, m1 in loader:
            t0 = t0.to(device)
            t1 = t1.to(device)
            m0 = m0.to(device)
            m1 = m1.to(device)
            offset0, offset1 = model.offsets()
            e0, w0, n0 = model.estimate_one(t0, m0, offset0)
            e1, w1, n1 = model.estimate_one(t1, m1, offset1)
            delta_timing.append((e0 - e1).cpu().numpy())

            for det, t, mask, weight, offset, neff in [
                (0, t0, m0, w0, offset0, n0),
                (1, t1, m1, w1, offset1, n1),
            ]:
                tc = t - offset
                big = torch.full_like(tc, 1.0e9)
                first_index = torch.argmin(torch.where(mask, tc, big), dim=1)
                first = torch.min(torch.where(mask, tc, big), dim=1).values
                dt = tc - first[:, None]

                weight_np = weight.cpu().numpy()
                dt_np = dt.cpu().numpy()
                first_np = first_index.cpu().numpy()
                mask_np = mask.cpu().numpy()

                weight_sum[det] += weight_np.sum(axis=0)
                neff_values[det].append(neff.cpu().numpy())

                for first_eoch, weights in zip(first_np, weight_np):
                    first_doch = int(EO2DO[first_eoch])
                    first_count[det, first_doch] += 1
                    for eoch, value in enumerate(weights):
                        weight_by_first_sum[det, first_doch, int(EO2DO[eoch])] += value

                bins = np.digitize(dt_np[mask_np], dt_edges) - 1
                values = weight_np[mask_np]
                valid = (bins >= 0) & (bins < nbins)
                np.add.at(dt_weight_sum[det], bins[valid], values[valid])
                np.add.at(dt_weight_count[det], bins[valid], 1)

    neff = np.stack([
        np.concatenate(neff_values[0]),
        np.concatenate(neff_values[1]),
    ])
    delta_timing = np.concatenate(delta_timing)

    avg_weight_eo = weight_sum / np.maximum(1, neff.shape[1])
    avg_weight_do = np.zeros_like(avg_weight_eo)
    for eoch, doch in enumerate(EO2DO):
        avg_weight_do[:, doch] = avg_weight_eo[:, eoch]

    avg_weight_by_first = np.divide(
        weight_by_first_sum,
        np.maximum(1, first_count)[:, :, None],
    )
    avg_weight_vs_dt = np.divide(
        dt_weight_sum,
        np.maximum(1, dt_weight_count),
    )

    return {
        "eo2do": EO2DO,
        "avg_weight_eo": avg_weight_eo,
        "avg_weight_do": avg_weight_do,
        "avg_weight_by_first_do": avg_weight_by_first,
        "first_count_do": first_count,
        "dt_edges": dt_edges,
        "avg_weight_vs_dt": avg_weight_vs_dt,
        "dt_weight_count": dt_weight_count,
        "neff": neff,
        "delta_timing": delta_timing,
    }


def write_diagnostic_plots(path, diagnostics):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    avg_weight_do = diagnostics["avg_weight_do"]
    weight_by_first = diagnostics["avg_weight_by_first_do"]
    first_count = diagnostics["first_count_do"]
    dt_edges = diagnostics["dt_edges"]
    avg_weight_vs_dt = diagnostics["avg_weight_vs_dt"]
    neff = diagnostics["neff"]
    delta_timing = diagnostics["delta_timing"]
    dt_center = 0.5 * (dt_edges[:-1] + dt_edges[1:])

    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
        for det in range(2):
            ax = axes[det, 0]
            image = avg_weight_do[det].reshape(8, 4)
            im = ax.imshow(image, origin="lower", aspect="auto")
            ax.set_title(f"TIMING{det}: average weight vs DO channel")
            ax.set_xlabel("DO x")
            ax.set_ylabel("DO y")
            fig.colorbar(im, ax=ax)

            ax = axes[det, 1]
            image = first_count[det].reshape(8, 4)
            im = ax.imshow(image, origin="lower", aspect="auto")
            ax.set_title(f"TIMING{det}: first-hit DO occupancy")
            ax.set_xlabel("DO x")
            ax.set_ylabel("DO y")
            fig.colorbar(im, ax=ax)
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        for det in range(2):
            im = axes[det].imshow(weight_by_first[det], origin="lower", aspect="auto")
            axes[det].set_title(f"TIMING{det}: avg weight DO vs first-hit DO")
            axes[det].set_xlabel("weighted DO channel")
            axes[det].set_ylabel("first-hit DO channel")
            fig.colorbar(im, ax=axes[det])
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        for det in range(2):
            axes[det].plot(dt_center, avg_weight_vs_dt[det], drawstyle="steps-mid")
            axes[det].set_title(f"TIMING{det}: average weight vs dt from first")
            axes[det].set_xlabel("dt from first calibrated hit [clock]")
            axes[det].set_ylabel("average assigned weight")
            axes[det].grid(True, alpha=0.3)
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        for det in range(2):
            axes[det].hist(neff[det], bins=60, range=(0, 32), histtype="step")
            axes[det].set_title(f"TIMING{det}: effective number of channels")
            axes[det].set_xlabel("1 / sum weights^2")
            axes[det].set_ylabel("events")
            axes[det].grid(True, alpha=0.3)
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        ax.hist(delta_timing, bins=2048, range=(-1.0, 1.0), histtype="step")
        ax.set_title("ML estimator delta timing")
        ax.set_xlabel("TIMING0 - TIMING1 [clock]")
        ax.set_ylabel("events")
        ax.grid(True, alpha=0.3)
        text = (
            f"entries = {len(delta_timing)}\n"
            f"mean = {np.mean(delta_timing):.6g}\n"
            f"RMS = {np.std(delta_timing):.6g}"
        )
        ax.text(0.98, 0.95, text, transform=ax.transAxes,
                ha="right", va="top")
        pdf.savefig(fig)
        plt.close(fig)


def format_stats(prefix, stats):
    return (
        f"{prefix} entries={stats['entries']} "
        f"rms={stats['rms']:.6g} mean={stats['mean']:.6g} "
        f"median={stats['median']:.6g} "
        f"q16={stats['q16']:.6g} q84={stats['q84']:.6g} "
        f"neff0={stats['neff0']:.3g} neff1={stats['neff1']:.3g}"
    )


def write_calibration(path, offset0, offset1):
    with open(path, "w", encoding="utf-8") as out:
        out.write("[CHANNEL]\n")
        out.write("# device fifo column pixel offset\n")
        out.write("# TIMING offsets learned by timing_ml.py\n")
        out.write("# calibrated_time = raw_time - offset\n")
        for det, offsets in enumerate([offset0, offset1]):
            fifo0 = 0 if det == 0 else 4
            fifo1 = 3 if det == 0 else 7
            for eoch, offset in enumerate(offsets):
                column = eoch // 4
                pixel = eoch % 4
                for fifo in range(fifo0, fifo1 + 1):
                    out.write(f"200 {fifo} {column} {pixel} {offset:.12g}\n")


def load_channel_offsets(path):
    offset = np.zeros((2, 32), dtype=np.float32)
    loaded = np.zeros((2, 32), dtype=bool)
    if not path:
        return offset, loaded

    with open(path, "r", encoding="utf-8") as inp:
        for line in inp:
            line = line.split("#", 1)[0].strip()
            if not line or line.startswith("["):
                continue
            parts = line.split()
            if len(parts) != 5:
                continue
            device, fifo, column, pixel = map(int, parts[:4])
            value = float(parts[4])
            if device != 200:
                continue
            det = 0 if 0 <= fifo <= 3 else 1 if 4 <= fifo <= 7 else -1
            if det < 0:
                continue
            eoch = pixel + 4 * column
            if 0 <= eoch < 32:
                offset[det, eoch] = value
                loaded[det, eoch] = True
    return offset, loaded


def main():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-root")
    source.add_argument("--input-npz")
    parser.add_argument("--output", default="timing_ml.npz")
    parser.add_argument("--model-output", default="")
    parser.add_argument("--load-model", default="")
    parser.add_argument("--diagnostics-output", default="")
    parser.add_argument("--diagnostics-plot-output", default="")
    parser.add_argument("--diagnostics-dt-range", type=float, default=2.0)
    parser.add_argument("--diagnostics-dt-bins", type=int, default=160)
    parser.add_argument("--calibration-output", default="timing_ml_offsets.conf")
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--lr", type=float, default=2.0e-3)
    parser.add_argument("--hidden", type=int, default=32)
    parser.add_argument("--temperature", type=float, default=0.12)
    parser.add_argument("--initial-calibration", default="")
    parser.add_argument("--neff-target", type=float, default=6.0)
    parser.add_argument("--neff-weight", type=float, default=1.0e-3)
    parser.add_argument("--offset-weight", type=float, default=1.0e-4)
    parser.add_argument("--max-events", type=int, default=0)
    parser.add_argument("--validation-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    if args.validation_fraction < 0.0 or args.validation_fraction >= 1.0:
        raise ValueError("--validation-fraction must satisfy 0 <= value < 1")
    if not args.model_output:
        args.model_output = str(Path(args.output).with_suffix(".pt"))

    if args.input_npz:
        t0, m0, t1, m1 = load_npz(args.input_npz)
    else:
        t0, m0, t1, m1 = load_root_with_uproot(args.input_root)

    good = m0.any(axis=1) & m1.any(axis=1)
    t0, m0, t1, m1 = t0[good], m0[good], t1[good], m1[good]
    if args.max_events > 0:
        t0, m0, t1, m1 = t0[: args.max_events], m0[: args.max_events], t1[: args.max_events], m1[: args.max_events]
    t0, t1 = normalize_event_times(t0, m0, t1, m1)

    rng = np.random.default_rng(args.seed)
    order = rng.permutation(len(t0))
    nvalid = int(round(args.validation_fraction * len(order)))
    valid_idx = order[:nvalid]
    train_idx = order[nvalid:]

    def make_dataset(indices):
        return TensorDataset(
            torch.tensor(t0[indices], dtype=torch.float32),
            torch.tensor(m0[indices], dtype=torch.bool),
            torch.tensor(t1[indices], dtype=torch.float32),
            torch.tensor(m1[indices], dtype=torch.bool),
        )

    train_dataset = make_dataset(train_idx)
    valid_dataset = make_dataset(valid_idx) if nvalid > 0 else train_dataset
    loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, drop_last=False)
    train_eval_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=False, drop_last=False)
    valid_eval_loader = DataLoader(valid_dataset, batch_size=args.batch_size, shuffle=False, drop_last=False)

    device = torch.device(args.device)
    checkpoint = None
    if args.load_model:
        checkpoint = torch.load(args.load_model, map_location=device, weights_only=False)
        args.hidden = int(checkpoint["hidden"])
        args.temperature = float(checkpoint["temperature"])

    initial, loaded = load_channel_offsets(args.initial_calibration)
    if checkpoint is not None:
        initial[0] = checkpoint["initial_offset0"].detach().cpu().numpy()
        initial[1] = checkpoint["initial_offset1"].detach().cpu().numpy()
        loaded[:, :] = True

    if args.initial_calibration and checkpoint is None:
        missing = np.argwhere(~loaded)
        if len(missing) > 0:
            print(f"WARNING: initial calibration is missing {len(missing)} TIMING entries")

    model = TimingEstimator(
        hidden=args.hidden,
        temperature=args.temperature,
        initial_offset0=torch.tensor(initial[0]),
        initial_offset1=torch.tensor(initial[1]),
    ).to(device)

    if checkpoint is not None:
        model.load_state_dict(checkpoint["state_dict"])

    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)

    for epoch in range(1, args.epochs + 1):
        model.train()
        total = 0.0
        nbatch = 0
        for bt0, bm0, bt1, bm1 in loader:
            bt0 = bt0.to(device)
            bt1 = bt1.to(device)
            bm0 = bm0.to(device)
            bm1 = bm1.to(device)
            e0, e1, _, _, neff0, neff1 = model(bt0, bm0, bt1, bm1)
            delta = e0 - e1
            delta_offset0, delta_offset1 = model.offset_deltas()
            loss = robust_loss(delta)
            loss = loss + args.neff_weight * (
                torch.relu(args.neff_target - neff0).mean()
                + torch.relu(args.neff_target - neff1).mean()
            )
            loss = loss + args.offset_weight * (delta_offset0.square().mean() + delta_offset1.square().mean())
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss.detach())
            nbatch += 1
        if epoch == 1 or epoch % 10 == 0 or epoch == args.epochs:
            train_stats = evaluate(model, train_eval_loader, device)
            valid_stats = evaluate(model, valid_eval_loader, device)
            print(f"epoch {epoch:4d} loss={total / max(1, nbatch):.6g}")
            print("  " + format_stats("train", train_stats))
            print("  " + format_stats("valid", valid_stats))

    train_stats = evaluate(model, train_eval_loader, device)
    valid_stats = evaluate(model, valid_eval_loader, device)
    offset0, offset1 = model.offsets()
    offset0 = offset0.detach().cpu().numpy()
    offset1 = offset1.detach().cpu().numpy()

    np.savez(
        args.output,
        offset0=offset0,
        offset1=offset1,
        train_stats=np.array([train_stats["entries"], train_stats["mean"], train_stats["rms"], train_stats["median"], train_stats["q16"], train_stats["q84"], train_stats["neff0"], train_stats["neff1"]]),
        valid_stats=np.array([valid_stats["entries"], valid_stats["mean"], valid_stats["rms"], valid_stats["median"], valid_stats["q16"], valid_stats["q84"], valid_stats["neff0"], valid_stats["neff1"]]),
    )
    torch.save(
        {
            "state_dict": model.state_dict(),
            "hidden": args.hidden,
            "temperature": args.temperature,
            "initial_offset0": model.initial_offset0.detach().cpu(),
            "initial_offset1": model.initial_offset1.detach().cpu(),
            "offset0": torch.tensor(offset0),
            "offset1": torch.tensor(offset1),
            "train_stats": train_stats,
            "valid_stats": valid_stats,
        },
        args.model_output,
    )
    write_calibration(args.calibration_output, offset0, offset1)

    if args.diagnostics_output:
        diagnostics = collect_diagnostics(
            model,
            valid_eval_loader,
            device,
            nbins=args.diagnostics_dt_bins,
            dt_range=(-args.diagnostics_dt_range, args.diagnostics_dt_range),
        )
        np.savez(args.diagnostics_output, **diagnostics)
        print(f"wrote {args.diagnostics_output}")
        if args.diagnostics_plot_output:
            write_diagnostic_plots(args.diagnostics_plot_output, diagnostics)
            print(f"wrote {args.diagnostics_plot_output}")

    print("final")
    print("  " + format_stats("train", train_stats))
    print("  " + format_stats("valid", valid_stats))
    print(f"wrote {args.output}")
    print(f"wrote {args.model_output}")
    print(f"wrote {args.calibration_output}")


if __name__ == "__main__":
    main()
