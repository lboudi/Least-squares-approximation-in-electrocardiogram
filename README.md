
# Least Squares Approximations for Electrocardiology

## Summary of the content

- **Objective**: To study how to numerically differentiate ECG signals, despite the noise and instability inherent in direct numerical differentiation.
- **Proposed method**: Use **polynomial approximations via the least squares method** to smooth the signals before computing their derivatives.
- **Technical approach**: From data points \((t_i, f_i)\), a polynomial \(p\) of fixed degree \(m\) is sought that minimizes the sum of squared errors:
  \[
  J(q) = \sum_i |q(x_i) - f_i|^2
  \]
  for all polynomials \(q\) of degree at most \(m\).
- **Application**: Study on real ECG signals, specifically in the case of atrial fibrillation (illustrated by a figure in the document).
