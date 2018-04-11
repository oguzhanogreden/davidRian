# News

## v0.9.1

- Approximate values in M matrix are replaced with analytically obtained values ([#9](https://github.com/oguzhanogreden/dcurver/issues/9)). Relevant code is reorganized.
- Gradient function raises a warning if `any(phi==0)` since the gradient returns 0 at this value.