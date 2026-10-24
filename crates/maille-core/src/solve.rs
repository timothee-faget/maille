use crate::assembly::DenseMatrix;
use crate::error::MailleError;

/// Résout A x = b par Cholesky (A symétrique définie positive).
/// Renvoie `Singular` si un pivot est négligeable devant la plus grande diagonale.
pub fn cholesky_solve(a: &DenseMatrix, b: &[f64]) -> Result<Vec<f64>, MailleError> {
    let n = a.n;
    let max_diag = (0..n).map(|i| a.get(i, i)).fold(0.0_f64, f64::max);
    let tol = 1e-12 * max_diag;

    // A = L Lᵀ, L triangulaire inférieure
    let mut l = DenseMatrix::zeros(n);
    for i in 0..n {
        for j in 0..=i {
            let s: f64 = (0..j).map(|k| l.get(i, k) * l.get(j, k)).sum();
            if i == j {
                let d = a.get(i, i) - s;
                if d.is_nan() || d <= tol {
                    return Err(MailleError::Singular);
                }
                l.add(i, i, d.sqrt());
            } else {
                l.add(i, j, (a.get(i, j) - s) / l.get(j, j));
            }
        }
    }

    // L y = b
    let mut x = vec![0.0; n];
    for i in 0..n {
        let s: f64 = (0..i).map(|k| l.get(i, k) * x[k]).sum();
        x[i] = (b[i] - s) / l.get(i, i);
    }
    // Lᵀ x = y
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|k| l.get(k, i) * x[k]).sum();
        x[i] = (x[i] - s) / l.get(i, i);
    }
    Ok(x)
}
