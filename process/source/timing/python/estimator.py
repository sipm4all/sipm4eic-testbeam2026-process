
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Sequence, Union
import json

import joblib
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

NATIVE_UNIT_NS = 3.125
NATIVE_UNIT_PS = 3125.0


class SigmaNet(nn.Module):
    def __init__(self, nin: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(nin, 64), nn.SiLU(),
            nn.Linear(64, 32), nn.SiLU(),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)


@dataclass
class TimingEstimate:
    T0_native: float
    sigma0_native: float
    T1_native: float
    sigma1_native: float
    T_event_native: float
    sigma_event_native: float

    @property
    def T0_ns(self): return self.T0_native * NATIVE_UNIT_NS
    @property
    def sigma0_ps(self): return self.sigma0_native * NATIVE_UNIT_PS
    @property
    def T1_ns(self): return self.T1_native * NATIVE_UNIT_NS
    @property
    def sigma1_ps(self): return self.sigma1_native * NATIVE_UNIT_PS
    @property
    def T_event_ns(self): return self.T_event_native * NATIVE_UNIT_NS
    @property
    def sigma_event_ps(self): return self.sigma_event_native * NATIVE_UNIT_PS

    def as_dict(self):
        d = asdict(self)
        d.update({
            "T0_ns": self.T0_ns,
            "sigma0_ps": self.sigma0_ps,
            "T1_ns": self.T1_ns,
            "sigma1_ps": self.sigma1_ps,
            "T_event_ns": self.T_event_ns,
            "sigma_event_ps": self.sigma_event_ps,
        })
        return d


class TimingEventEstimator:
    """
    Final trained estimator.

    Input:
      columns 0..31  = TIMING0 DO0..DO31
      columns 32..63 = TIMING1 DO0..DO31

    Input times must already be calibrated with:
      calibrated_time = hit.time - offset_860
    """

    def __init__(self, model_dir: Union[str, Path] = None, device: str = "cpu"):
        self.model_dir = Path(model_dir) if model_dir else Path(__file__).resolve().parent
        self.device = torch.device(device)

        # Model B / sigma feature scalers.
        self.sc0 = joblib.load(self.model_dir / "scaler_timing0.joblib")
        self.sc1 = joblib.load(self.model_dir / "scaler_timing1.joblib")

        # Two independently trained time estimators that are blended.
        self.timeA0 = joblib.load(self.model_dir / "timeA_timing0.joblib")
        self.timeA1 = joblib.load(self.model_dir / "timeA_timing1.joblib")
        self.timeB0 = joblib.load(self.model_dir / "timeB_timing0.joblib")
        self.timeB1 = joblib.load(self.model_dir / "timeB_timing1.joblib")

        nin = int(self.sc0.n_features_in_)
        self.snet0 = SigmaNet(nin).to(self.device)
        self.snet1 = SigmaNet(nin).to(self.device)
        self.snet0.load_state_dict(torch.load(self.model_dir / "sigmanet_timing0.pt", map_location=self.device))
        self.snet1.load_state_dict(torch.load(self.model_dir / "sigmanet_timing1.pt", map_location=self.device))
        self.snet0.eval()
        self.snet1.eval()

        with open(self.model_dir / "metadata.json") as f:
            self.meta = json.load(f)
        self.a = float(self.meta["time_model"]["ensemble_blend_model_A"])
        self.sigma_scale = float(self.meta["sigma_model"]["variance_calibration_scale"])

    @staticmethod
    def _anchor_and_sorted(t):
        t = np.asarray(t, dtype=np.float64)
        s = np.sort(t, axis=1)
        B = s[:, :10].mean(axis=1)
        return B, s

    @staticmethod
    def _features_A(t, s, B):
        # Exact feature set used by final Model A: 62 features.
        rel = t - B[:, None]                 # 32
        gaps = s[:, :24] - B[:, None]       # 24
        summaries = np.column_stack([
            s[:,0] - B,
            s[:,1] - s[:,0],
            s[:,3] - s[:,0],
            s[:,7] - s[:,0],
            s[:,15] - s[:,0],
            s[:,23] - s[:,0],
        ])                                   # 6
        return np.column_stack([rel, gaps, summaries]).astype(np.float32)

    @staticmethod
    def _features_B(t, s, B):
        # Exact feature set used by Model B and sigma networks: 52 features.
        rel = t - B[:, None]                 # 32
        gaps = s[:, :20] - B[:, None]       # 20
        return np.column_stack([rel, gaps]).astype(np.float32)

    def _predict_detector(self, t, detector):
        t = np.asarray(t, dtype=np.float64)
        B, s = self._anchor_and_sorted(t)

        XA = self._features_A(t, s, B)
        XB = self._features_B(t, s, B)

        if detector == 0:
            XBs = self.sc0.transform(XB).astype(np.float32)
            pA = self.timeA0.predict(XA)
            pB = self.timeB0.predict(XBs)
            T = B - (self.a*pA + (1.0-self.a)*pB)
            snet = self.snet0
        else:
            XBs = self.sc1.transform(XB).astype(np.float32)
            pA = self.timeA1.predict(XA)
            pB = self.timeB1.predict(XBs)
            T = B + (self.a*pA + (1.0-self.a)*pB)
            snet = self.snet1

        xt = torch.from_numpy(XBs).to(self.device)
        with torch.no_grad():
            raw_sigma = torch.sqrt(F.softplus(snet(xt)) + 1e-6).cpu().numpy()
        sigma = raw_sigma * self.sigma_scale
        return T.astype(np.float64), sigma.astype(np.float64)

    def predict_batch(self, events):
        arr = np.asarray(events, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.ndim != 2 or arr.shape[1] != 64:
            raise ValueError(f"Expected shape (N,64), got {arr.shape}")
        if not np.isfinite(arr).all():
            raise ValueError("Input contains NaN or inf.")

        T0, s0 = self._predict_detector(arr[:, :32], 0)
        T1, s1 = self._predict_detector(arr[:, 32:], 1)

        Te = 0.5*(T0+T1)
        se = 0.5*np.sqrt(s0*s0+s1*s1)

        return {
            "T0_native":T0, "sigma0_native":s0,
            "T1_native":T1, "sigma1_native":s1,
            "T_event_native":Te, "sigma_event_native":se,
            "T0_ns":T0*NATIVE_UNIT_NS, "sigma0_ps":s0*NATIVE_UNIT_PS,
            "T1_ns":T1*NATIVE_UNIT_NS, "sigma1_ps":s1*NATIVE_UNIT_PS,
            "T_event_ns":Te*NATIVE_UNIT_NS, "sigma_event_ps":se*NATIVE_UNIT_PS,
        }

    def predict(self, event: Sequence[float]) -> TimingEstimate:
        o = self.predict_batch(event)
        return TimingEstimate(
            T0_native=float(o["T0_native"][0]),
            sigma0_native=float(o["sigma0_native"][0]),
            T1_native=float(o["T1_native"][0]),
            sigma1_native=float(o["sigma1_native"][0]),
            T_event_native=float(o["T_event_native"][0]),
            sigma_event_native=float(o["sigma_event_native"][0]),
        )
