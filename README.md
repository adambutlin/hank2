# hank2 — Wealth and Consumption in the UK: A Two-Asset Heterogeneous Agent Model

MATLAB code and figures accompanying the paper:

> **Butlin, A. (2024). "Wealth and Consumption in the UK: A Two-Asset Heterogeneous Agent Model."**
> *UCL Journal of Economics*, vol. 3 no. 1. DOI: [10.14324/111.444.2755-0877.1889](https://doi.org/10.14324/111.444.2755-0877.1889)

The full paper is included as [`uje-1889-butlin.pdf`](uje-1889-butlin.pdf).

## Overview

This project calibrates a **two-asset heterogeneous agent model** with uninsurable
idiosyncratic income risk and no aggregate uncertainty to the UK income and wealth
distribution, using household balance-sheet data from the Office for National
Statistics' Wealth and Assets Survey (WAS).

Households hold two assets:

- a **liquid** asset with a lower return, and
- an **illiquid** asset with a higher return whose withdrawals incur a transaction cost.

This creates a trade-off between holding liquid wealth to self-insure against income
shocks and earning the higher illiquid return. The model extends the classic
Bewley–Huggett–Aiyagari framework and is used to study how return differentials and
transaction costs shape the stationary distributions of income and wealth, and what
this implies for the effectiveness of fiscal and monetary policy in the UK.

## Repository contents

### Model code

| File | Description |
|------|-------------|
| [`two_asset.m`](two_asset.m) | Main solver for the two-asset (liquid/illiquid) heterogeneous agent model |
| [`aiyagari.m`](aiyagari.m) | One-asset Aiyagari model (benchmark) |
| [`two_asset_kinked_cost.m`](two_asset_kinked_cost.m) | Kinked illiquid-asset adjustment cost function |
| [`two_asset_kinked_FOC.m`](two_asset_kinked_FOC.m) | First-order condition for the kinked adjustment cost |
| [`test.m`](test.m) | Scratch / testing script |

### Figures

The `.eps` files are the figures from the paper:

| File | Paper figure |
|------|--------------|
| [`consumption.eps`](consumption.eps) | Figure 1 — Consumption policy for low and high types |
| [`deposits.eps`](deposits.eps) | Figure 2 — Deposits (illiquid contributions) for low and high types |
| [`sav.eps`](sav.eps) | Figure 3 — Liquid savings for low and high types |
| [`ill_sav.eps`](ill_sav.eps) | Figure 4 — Illiquid savings for low and high types |
| [`stat_dist.eps`](stat_dist.eps) | Figure 5 — Stationary distribution for low and high types |
| [`adj_cost.eps`](adj_cost.eps) | Illiquid-asset adjustment cost schedule |

## Running

The scripts are written for MATLAB. Open and run [`two_asset.m`](two_asset.m) in MATLAB
to solve the two-asset model and reproduce the policy functions and stationary
distribution; run [`aiyagari.m`](aiyagari.m) for the one-asset benchmark.

## Author

Adam Butlin — MSc Economics (2022–23), Department of Economics, University College London.

## Licence

The paper is open access under the [Creative Commons Attribution Licence (CC BY) 4.0](https://creativecommons.org/licenses/by/4.0/).
