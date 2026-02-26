# -*- coding: utf-8 -*-
"""
SeqGen_Feasibility_space.py (Z-only loss)
----------------------------------------
Annealed swap optimizer that minimizes squared error in Z-space:
    E = 0.5*(Z_p - zP_target)^2 + 0.5*(Z_H - zH_target)^2
"""
import numpy as np
import random
from math import exp
from typing import List, Tuple
from PatternMetric import StickerPattern
from NormMetric import normalised_H, normalised_P

# Matplotlib only used in plot_seq; safe to import lazily there if desired.
import matplotlib.pyplot as plt

class SeqDesign:
    def __init__(self, n_sticker: int, length: int, htarget: float, ptarget: float, sticker_char: str = 'S'):
        self.n_sticker = int(n_sticker)
        self.length = int(length)
        self.H_target = float(htarget)   # interpreted as Z(H) target
        self.p_target = float(ptarget)   # interpreted as Z(p) target
        self.sticker = sticker_char

    # --------------------- Core objective in Z-space ---------------------
    def energy(self, seq: List[str]) -> float:
        s = ''.join(seq)
        z_p = normalised_P(s, sticker=self.sticker)
        z_h = normalised_H(s, sticker=self.sticker)
        return 0.5 * (z_p - self.p_target) ** 2 + 0.5 * (z_h - self.H_target) ** 2

    # ----------------------- Initial sequence ---------------------------
    def _init_sequence(self) -> List[str]:
        nS = self.n_sticker
        nX = self.length - nS
        if nS < 0 or nX < 0:
            raise ValueError("Invalid (nsticker, length). nsticker must be in [0, length].")
        seq = ['S'] * nS + ['X'] * nX
        random.shuffle(seq)
        return seq

    # ----------------------- Single annealing run -----------------------
    def run_monte_carlo(self, seq, steps: int, T_init: float = 10.0, T_decay: float = 0.999,
                    log_steps: bool = False, step_log_frac: float = 0.20):
        """
        If log_steps is True, print a one-line update every ~step_log_frac of steps.
        """
        import random
        from math import exp

        current_seq = seq[:]
        current_cost = self.energy(current_seq)
        best_seq = current_seq[:]
        best_cost = current_cost
        T = T_init

        # Ensure each step proposes a real swap: pick from S and X pools
        S_idx = [i for i, c in enumerate(current_seq) if c == self.sticker]
        X_idx = [i for i, c in enumerate(current_seq) if c != self.sticker]
        if len(S_idx) == 0 or len(X_idx) == 0:
            return best_seq, best_cost

        # Step logging cadence
        k = None
        if log_steps and steps >= 5:
            k = max(1, int(steps * float(step_log_frac)))

        for step in range(1, steps + 1):
            i = random.choice(S_idx)
            j = random.choice(X_idx)
            if i == j:
                T *= T_decay
                continue

            # Propose swap
            current_seq[i], current_seq[j] = current_seq[j], current_seq[i]
            new_cost = self.energy(current_seq)
            dE = new_cost - current_cost

            if dE <= 0 or random.random() < exp(-dE / max(T, 1e-12)):
                current_cost = new_cost
                # update pools after swap
                if current_seq[i] == self.sticker:
                    if i not in S_idx: S_idx.append(i)
                    if j not in X_idx: X_idx.append(j)
                    if i in X_idx: X_idx.remove(i)
                    if j in S_idx: S_idx.remove(j)
                else:
                    if i not in X_idx: X_idx.append(i)
                    if j not in S_idx: S_idx.append(j)
                    if i in S_idx: S_idx.remove(i)
                    if j in X_idx: X_idx.remove(j)
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_seq = current_seq[:]
            else:
                # reject → revert
                current_seq[i], current_seq[j] = current_seq[j], current_seq[i]

            T *= T_decay

            if k and (step % k == 0 or step == steps):
                pct = int(round(100 * step / steps))
                print(f"    step {step}/{steps} ({pct:3d}%) | current={current_cost:.6g} | best={best_cost:.6g}",
                    flush=True)

        return best_seq, best_cost


    # ------------------------- Multi-epoch driver -----------------------
    def generate_seq(self, nStep: int = 25, steps: int = 50000,
                 log_epochs: bool = True, log_steps: bool = False, step_log_frac: float = 0.20):
        import numpy as np, time
        seq_init = self._init_sequence()
        costs = []
        best_overall = seq_init[:]
        best_cost = self.energy(best_overall)

        t0 = time.time()
        for epoch in range(1, nStep + 1):
            best_seq_epoch, energy_epoch = self.run_monte_carlo(
                seq_init, steps=steps, T_init=10.0, T_decay=0.999,
                log_steps=log_steps, step_log_frac=step_log_frac
            )
            seq_init = best_seq_epoch
            costs.append(energy_epoch)
            if energy_epoch < best_cost:
                best_cost = energy_epoch
                best_overall = best_seq_epoch[:]

            if log_epochs:
                elapsed = time.time() - t0
                print(f"[{epoch}/{nStep}] epochs complete | best_cost={best_cost:.6g} "
                    f"| epoch_cost={energy_epoch:.6g} | elapsed={elapsed:.1f}s",
                    flush=True)

        step_tot = np.arange(steps, nStep * steps + 1, steps)
        return ''.join(best_overall), step_tot, costs


    # --------------------------- Diagnostics ----------------------------
    def plot_seq(self, best_seq: str, steps=None, E=None):
        # Compute Z-metrics for the best sequence
        z_p = normalised_P(best_seq, sticker=self.sticker)
        z_h = normalised_H(best_seq, sticker=self.sticker)

        fig, axs = plt.subplots(1, 2, figsize=(8, 4))
        # Left: sequence bead plot
        axs[0].set_title(f"Best seq (Zp={z_p:.2f}, Zh={z_h:.2f})")
        xs = np.arange(len(best_seq))
        ys = np.zeros_like(xs, dtype=float)
        colors = ['tab:blue' if c == 'S' else 'tab:gray' for c in best_seq]
        axs[0].scatter(xs, ys, s=20, c=colors)
        axs[0].set_yticks([])
        axs[0].set_xlabel('Position')

        # Right: cost curve
        if steps is not None and E is not None:
            axs[1].plot(steps, E, lw=1.5)
            axs[1].set_xlabel('Steps')
            axs[1].set_ylabel('Energy (Z-space)')
            axs[1].set_title('Convergence')
        plt.tight_layout()
        plt.show()
