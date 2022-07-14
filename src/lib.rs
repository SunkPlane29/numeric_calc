mod file;
mod root;
pub use file::DataFile;
pub use root::*;

// calculates the derivative of a certain function at a certain point, with a
// number of steps that defaults to 20, this doesn't affect much in terms of performance,
// but there is no need to have like 100 steps
pub fn calculate_derivative(fx: fn(f64) -> f64, x: f64) -> f64 {
    let eps = 1.0/(2.0_f64).powf(20.0); // I think this is a magic number for me hehe

    (fx(x + eps) - fx(x - eps)) / (2.0 * eps)
}

// calculates the riemann sum of a certain function at a certain point with an optinal number of
// rectangles which defaults to 100000, this is a very slow function, but it is a good way to
// calculate the area of a function
pub fn calculate_integral_with_rectangles(fx: fn(f64) -> f64, xmin: f64, xmax: f64, n_rectangles: Option<i64>) -> f64 {
    let n_rectangles = n_rectangles.unwrap_or(100000);

    let delx = (xmax - xmin)/n_rectangles as f64;

    let mut integral = 0.0;

    for i in 0..n_rectangles {
        let xmin_temp = xmin + (i as f64 * delx);
        let xmax_temp = xmin + ((i+1) as f64 * delx);

        integral += (xmax_temp - xmin_temp)*fx(xmin_temp + (delx/2.0));
    }

    integral
}

// --------------------------------------------------------------------------------------------------------------------
//
// CRATE TESTS
//
// --------------------------------------------------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_calculate_derivative() {
        struct Test {
            fx: fn(f64) -> f64,
            x: f64,
            expect: f64,
        }

        let tests = vec![
            Test{
                fx: |x| {x*2.0},
                x: 45.0,
                expect: 2.0,
            },
            Test{
                fx: |x| {x.powf(2.0)},
                x: 5.0,
                expect: 10.0,
            },
            Test{
                fx: |x| {x.powf(3.0)},
                x: 2.0,
                expect: 12.0,
            },
            Test{
                fx: |x| {x.powf(4.0)},
                x: 1.0,
                expect: 4.0,
            },
        ];

        let precision = 1.0e-16;

        for test in tests {
            let derivative = calculate_derivative(test.fx, test.x);
            assert!((derivative - test.expect).abs() < precision);
        }
    }

    #[test]
    fn test_calculate_integral_with_rectangles() {
        struct Test {
            fx: fn(f64) -> f64,
            xmin: f64,
            xmax: f64,
            expect: f64,
        }

        let tests = vec![
            Test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xmin: 0.0,
                xmax: 4.0,
                expect: 5.333333,
            },
            Test{
                fx: |x| {x.powf(2.0) + x - 6.0},
                xmin: -5.0,
                xmax: -1.0,
                expect: 5.333333,
            },
            Test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 0.5,
                xmax: 1.5,
                expect: -0.25,
            },
            Test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 1.5,
                xmax: 2.5,
                expect: 1.00614e-16,
            },
            Test{
                fx: |x| {(x.powf(3.0)) + (-6.0*x.powf(2.0)) + (11.0*x) - 6.0},
                xmin: 2.5,
                xmax: 3.5,
                expect: 0.25,
            },
        ];

        let n_rectangles = 100000;
        let precision = 1.0e-6;

        for test in tests {
            let integral = calculate_integral_with_rectangles(test.fx, test.xmin, test.xmax, Some(n_rectangles));
            assert!((integral - test.expect).abs() < precision);
        }
    }
}
