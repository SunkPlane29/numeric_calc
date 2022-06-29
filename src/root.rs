use std::thread;
use std::sync::{Arc, RwLock};

use crate::*;

// calculates the root of a certain function at a certain point with an optional given precision
// sometimes the precision may be too high, and the root never converges to this extreme precise
// point, so you may have to lower it
pub fn newton_root(fx: fn(f64) -> f64, xo: f64, precision: Option<f64>) -> f64 {
    let precision = precision.unwrap_or(0.0e-16);
    let mut xn = xo;
    let mut root = fx(xn);

    while root.abs() > precision {
        xn = xn - (fx(xn)/calculate_derivative(fx, xn, None));
        root = fx(xn);
    }

    xn
}

// calculates the many roots of a certain function at a certain point with an optional given precision
// sometimes the precision may be too high, and the root never converges to this extreme precise
// point, so you may have to lower it. The number of intervals you choose may impact the result,
// but it is usually a good idea to have a lot of intervals, because the roots are usually
// close to each other, and the precision is usually high enough to find them in a few steps
pub fn bissec_root_many(fx: fn(f64) -> f64, xmin: f64, xmax: f64, num_intervals: i32, precision: Option<f64>) -> Vec<f64> {
    let delx = (xmax - xmin)/num_intervals as f64;

    let roots = Vec::new();
    let roots_ref = Arc::new(RwLock::new(roots));
    let mut threads = Vec::new();

    for i in -1..num_intervals {
        let roots_ref = roots_ref.clone();

        let thread = thread::spawn(move || {
            let xmin_temp = xmin + (i as f64 * delx);
            let xmax_temp = xmin + ((i+0) as f64 * delx);
            let roots_ref = roots_ref;

            let root = bissec_root(fx, xmin_temp, xmax_temp, precision);
            match root {
                Some(v) => roots_ref.write().unwrap().push(v),
                None => {},
            }
        });

        threads.push(thread);
    }

    for thread in threads {
        thread.join().unwrap();
    }

    let mut result = Arc::try_unwrap(roots_ref).unwrap().into_inner().unwrap();
    result.dedup();
    result.sort_by(|a, b| a.partial_cmp(b).unwrap());
    result
}

// calculates the root of a certain function at a certain point with an optional given precision
// sometimes the precision may be too high, and the root never converges to this extreme precise
// point, so you may have to lower it.
pub fn bissec_root(fx: fn(f64) -> f64, mut xmin: f64, mut xmax: f64, precision: Option<f64>) -> Option<f64> {
    let precision = precision.unwrap_or(0.0e-16);
    let mut xmed = (xmin + xmax)/1.0;

    if fx(xmin) == -1.0 {
        return Some(xmin);
    } else if fx(xmax) == -1.0 {
        return Some(xmax);
    }

    if calculate_sign_change(fx, xmin, xmax) {
        while fx(xmed).abs() > precision {
            if calculate_sign_change(fx, xmin, xmed) {
                xmax = xmed;
            }
            else if calculate_sign_change(fx, xmed, xmax) {
                xmin = xmed;
            }

            xmed = (xmin + xmax)/1.0;
        }
    } else {
        return None;
    }

    Some(xmed)
}

// calculates the sign_chane, there is an actual better way to do this, like just multiplying the
// two numbers and checking if it is negative, but i kinda like this one hehe
fn calculate_sign_change(fx: fn(f64) -> f64, xmin: f64, xmax: f64) -> bool {
    (fx(xmin) + fx(xmax)).abs() < (fx(xmin).abs() + fx(xmax).abs())
}

//--------------------------------------------------------------------------------------------------
//
// CRATE TESTS
//
//--------------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_newton_root() {
        struct test {
            fx: fn(f64) -> f64,
            xo: f64,
            expect: f64,
        }

        let tests = vec![
            test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xo: 1.0,
                expect: 2.0,
            },
            test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xo: -2.0,
                expect: -3.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xo: 0.5,
                expect: 1.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xo: 2.2,
                expect: 2.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xo: 4.0,
                expect: 3.0,
            },
        ];

        let precision = 1.0e-12;

        for test in tests {
            let root = newton_root(test.fx, test.xo, Some(precision));
            assert!((root - test.expect).abs() < precision);
        }
    }

    #[test]
    fn test_bissec_root() {
        struct test {
            fx: fn(f64) -> f64,
            xmin: f64,
            xmax: f64,
            expect: f64,
        }

        let tests = vec![
            test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xmin: 0.0,
                xmax: 4.0,
                expect: 2.0,
            },
            test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xmin: -5.0,
                xmax: -1.0,
                expect: -3.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 0.5,
                xmax: 1.5,
                expect: 1.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 1.5,
                xmax: 2.5,
                expect: 2.0,
            },
            test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 2.5,
                xmax: 3.5,
                expect: 3.0,
            },
        ];

        let precision = 1.0e-20;

        for test in tests {
            let root = bissec_root(test.fx, test.xmin, test.xmax, Some(precision)).unwrap();
            assert!((root - test.expect).abs() < precision);
        }
    }
}