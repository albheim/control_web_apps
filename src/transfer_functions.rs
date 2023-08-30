#![allow(non_snake_case)]

pub trait TransferFunction {
    fn step_response(&self, t: f64) -> f64;
    fn bode_amplitude(&self, w: f64) -> f64;
    fn bode_phase(&self, w: f64) -> f64;
    fn poles(&self) -> Vec<[f64; 2]>;
    fn adjust_poles_to(&mut self, re: f64, im: f64);
}

#[derive(Debug, Clone, Copy)]
pub struct FirstOrderSystem {
    // first order system 1/(sT + 1)
    // pole = -1/T
    // https://www.tutorialspoint.com/control_systems/control_systems_response_first_order.htm
    pub T: f64,
    pub T_lower: f64,
    pub T_upper: f64,
}

impl TransferFunction for FirstOrderSystem {
    fn poles(&self) -> Vec<[f64; 2]> {
        vec![[-1.0 / self.T, 0.0]]
    }

    fn step_response(&self, t: f64) -> f64 {
        if t >= 0.0 {
            1.0 - (-t / self.T).exp()
        } else {
            0.0
        }
    }

    fn bode_amplitude(&self, w: f64) -> f64 {
        1.0 / (((w * self.T).powi(2) + 1.0).sqrt())
    }

    fn bode_phase(&self, w: f64) -> f64 {
        -(w * self.T).atan()
    }

    fn adjust_poles_to(&mut self, re: f64, _im: f64) {
        let pole_bound = -1.0/self.T_upper;

        if re >= pole_bound {
            self.T = -1.0 / pole_bound;
        } else {
            self.T = -1.0 / re;
        }
    }
}



#[derive(Debug, Clone, Copy)]
pub struct SecondOrderSystem {
    // second order system w^2/(s^2 + 2dw s + w^2)
    // poles = -dw +- w sqrt(d^2 - 1)
    // https://www.tutorialspoint.com/control_systems/control_systems_response_second_order.htm
    pub d: f64,
    pub w: f64,
    pub d_lower: f64,
    pub d_upper: f64,
    pub w_lower: f64,
    pub w_upper: f64,
}

impl TransferFunction for SecondOrderSystem {
    fn poles(&self) -> Vec<[f64; 2]> {
        let (d, w) = (self.d, self.w);

        if d == 0.0 {
            vec![[0.0, w], [0.0, -w]]
        } else if (0.0 < d) && (d < 1.0) {
            vec![
                [-d * w, (1.0 - d.powi(2)).sqrt() * w],
                [-d * w, -(1.0 - d.powi(2)).sqrt() * w],
            ]
        } else if d == 1.0 {
            vec![[-w, 0.0], [-w, 0.0]]
        } else {
            vec![
                [-d * w + (d.powi(2) - 1.0).sqrt() * w, 0.0],
                [-d * w - (d.powi(2) - 1.0).sqrt() * w, 0.0],
            ]
        }
    }

    fn step_response(&self, t: f64) -> f64 {
        let (d, w) = (self.d, self.w);

        if t < 0.0 {
            return 0.0;
        }

        if d == 0.0 {
            1.0 - (w * t).cos()
        } else if (0.0 < d) && (d < 1.0) {
            let d_1_sqrt = (1.0 - d.powi(2)).sqrt();
            let w_d = w*d_1_sqrt;
            1.0 - ( (-d*w*t).exp() ) * ( (w_d*t).cos() ) - ( (-d*w*t).exp() ) * ( (w_d*t).sin() ) * d / d_1_sqrt
        } else if d == 1.0 {
            1.0 - ((-w * t).exp()) * (1.0 + w * t)
        } else {
            let d_1_sqrt = (d.powi(2) - 1.0).sqrt();
            1.0 + (-t * w * (d + d_1_sqrt)).exp() / (2.0 * (d + d_1_sqrt) * d_1_sqrt)
                - (-t * w * (d - d_1_sqrt)).exp() / (2.0 * (d - d_1_sqrt) * d_1_sqrt)
        }
    }

    fn bode_amplitude(&self, w: f64) -> f64 {
        let (d, wp) = (self.d, self.w);

        wp.powi(2) / ( ( (wp.powi(2) - w.powi(2)).powi(2) + (2f64*d*wp*w).powi(2) ).sqrt() )
    }

    fn bode_phase(&self, w: f64) -> f64 {
        use std::f64::consts::PI;

        let (d, wp) = (self.d, self.w);

        let ph = -( (2f64*d*wp*w)/(wp.powi(2) - w.powi(2))  ).atan();
        if ph > 0.0 {
            ph - PI
        } else {
            ph
        }
    }

    fn adjust_poles_to(&mut self, re: f64, im: f64) {
        if re >= 0.0 {
            return
        }

        let mut d_new = self.d;
        let w_new;

        if self.d < 1.0 {
            // two complex poles
            let (re2, im2) = (re.powi(2), im.powi(2));
            d_new = (re2/(re2+im2)).sqrt();
            w_new = -re/d_new;
        } else if self.d == 1.0 {
            w_new = -re;
            // real double pole
        } else {
            // two real poles
            let d2 = self.d.powi(2);
            let fast = -self.d*self.w + self.w*(d2 - 1.0).sqrt();
            let slow = -self.d*self.w - self.w*(d2 - 1.0).sqrt();

            if (re - fast).abs() < (re - slow).abs() {
                w_new = re/( -self.d + (d2-1.0).sqrt() );
            } else {
                w_new = re/( -self.d - (d2-1.0).sqrt() );
            }
        }

        if self.d_lower <= d_new && self.d_upper >= d_new && self.w_lower <= w_new && self.w_upper >= w_new {
            self.d = d_new;
            self.w = w_new;
        }
    }
}
