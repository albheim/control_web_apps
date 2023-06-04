#![allow(non_snake_case)]

pub trait TransferFunction {
    fn step_response(&self, t: f64) -> f64;
    fn bode_amplitude(&self, w: f64) -> f64;
    fn bode_phase(&self, w: f64) -> f64;
    fn poles(&self) -> Vec<[f64; 2]>;
    fn zeros(&self) -> Vec<[f64; 2]> { vec![] }
    fn adjust_pole_zero(&mut self, re: f64, im: f64);
}

#[derive(Debug, Clone, Copy)]
pub struct FirstOrderSystem {
    // first order system exp(-sL)*K/(sT + 1)
    // pole = -1/T
    // https://www.tutorialspoint.com/control_systems/control_systems_response_first_order.htm
    pub T: f64,
    pub K: f64,
    pub L: f64,
    pub T_lower: f64,
    pub T_upper: f64,
    pub K_lower: f64,
    pub K_upper: f64,
    pub L_lower: f64,
    pub L_upper: f64,
}

impl TransferFunction for FirstOrderSystem {
    fn poles(&self) -> Vec<[f64; 2]> {
        vec![[-1.0 / self.T, 0.0]]
    }

    fn step_response(&self, t: f64) -> f64 {
        let t = t - self.L;
        self.K * if t >= 0.0 {
            1.0 - (-t / self.T).exp()
        } else {
            0.0
        }
    }

    fn bode_amplitude(&self, w: f64) -> f64 {
        self.K.abs() / (((w * self.T).powi(2) + 1.0).sqrt())
    }

    fn bode_phase(&self, w: f64) -> f64 {
        -self.L * w - (w * self.T).atan()
    }

    fn adjust_pole_zero(&mut self, re: f64, _im: f64) {
        let pole_bound = -1.0/self.T_upper;

        if re >= pole_bound {
            self.T = -1.0 / pole_bound;
        } else {
            self.T = -1.0 / re;
        }
    }
}



#[derive(Debug, Clone, Copy)]
pub struct RealSecondOrderSystem {
    // double first order system K * exp(-sL) * / ((sT1+1)*(sT2+1))
    // poles = -1/T_1, -1/T_2
    // https://www.tutorialspoint.com/control_systems/control_systems_response_first_order.htm
    pub T1: f64,
    pub T2: f64,
    pub K: f64,
    pub L: f64,
    pub T_lower: f64,
    pub T_upper: f64,
    pub K_lower: f64,
    pub K_upper: f64,
    pub L_lower: f64,
    pub L_upper: f64,
}

impl TransferFunction for RealSecondOrderSystem {
    fn poles(&self) -> Vec<[f64; 2]> {
        vec![[-1.0 / self.T1, 0.0], [-1.0 / self.T2, 0.0]]
    }

    fn step_response(&self, t: f64) -> f64 {
        let t = t - self.L;
        self.K * if t >= 0.0 {
            1.0 - (-self.T1 * (-t / self.T1).exp() + self.T2 * (-t / self.T2).exp()) / (self.T2 - self.T1)
        } else {
            0.0
        }
    }

    fn bode_amplitude(&self, w: f64) -> f64 {
        self.K.abs() / (((w * self.T1).powi(2) + 1.0).sqrt() * ((w * self.T2).powi(2) + 1.0).sqrt())
    }

    fn bode_phase(&self, w: f64) -> f64 {
        -self.L * w - (w * self.T1).atan() - (w * self.T2).atan()
    }

    fn adjust_pole_zero(&mut self, re: f64, _im: f64) {
        let pole_bound = -1.0/self.T_upper;
        let x1 = -1.0/self.T1;
        let x2 = -1.0/self.T2;

        let target = if re >= pole_bound {
            -1.0 / pole_bound
        } else {
            -1.0 / re
        };
        if (re - x1).abs() < (re - x2).abs() {
            self.T1 = target;
        } else {
            self.T2 = target;
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub struct RealSecondOrderSystemAndZero {
    // second order system exp(-sL) * K * (sT + 1) * w^2/(s^2 + 2dw s + w^2)
    // poles = -dw +- w sqrt(d^2 - 1)
    // https://www.tutorialspoint.com/control_systems/control_systems_response_second_order.htm
    pub T1: f64,
    pub T2: f64,
    pub Tz: f64,
    pub K: f64,
    pub L: f64,
    pub T_lower: f64,
    pub T_upper: f64,
    pub Tz_lower: f64,
    pub Tz_upper: f64,
    pub K_lower: f64,
    pub K_upper: f64,
    pub L_lower: f64,
    pub L_upper: f64,
}

impl TransferFunction for RealSecondOrderSystemAndZero {
    fn poles(&self) -> Vec<[f64; 2]> {
        vec![[-1.0 / self.T1, 0.0], [-1.0 / self.T2, 0.0]]
    }

    fn zeros(&self) -> Vec<[f64; 2]> {
        vec![[-1.0 / self.Tz, 0.0]]
    }

    fn step_response(&self, t: f64) -> f64 {
        let t = t - self.L;
        let (T1, T2, Tz) = (self.T1, self.T2, self.Tz);
        self.K * if t >= 0.0 {
            1.0 - ((Tz - T1) * (-t / T1).exp() + (T2 - Tz) * (-t / T2).exp()) / (T2 - T1)
        } else {
            0.0
        }
    }

    fn bode_amplitude(&self, w: f64) -> f64 {
        self.K.abs() * ((w * self.Tz).powi(2) + 1.0).sqrt() / (((w * self.T1).powi(2) + 1.0).sqrt() * ((w * self.T2).powi(2) + 1.0).sqrt())
    }

    fn bode_phase(&self, w: f64) -> f64 {
        -self.L * w + (w * self.Tz).atan() - (w * self.T1).atan() - (w * self.T2).atan()
    }

    fn adjust_pole_zero(&mut self, re: f64, _im: f64) {
        let x1 = -1.0/self.T1;
        let x2 = -1.0/self.T2;
        let xz = -1.0/self.Tz;

        if (re - xz).abs() < (re - x1).abs() && (re - xz).abs() < (re - x2).abs() {
            self.Tz = -1.0 / re
        } else if (re - x1).abs() < (re - x2).abs() {
            self.T1 = if self.T_upper * re >= -1.0 { self.T_upper } else { -1.0 / re }
        } else {
            self.T2 = if self.T_upper * re >= -1.0 { self.T_upper } else { -1.0 / re }
        }
    }
}




#[derive(Debug, Clone, Copy)]
pub struct ComplexSecondOrderSystem {
    // second order system exp(-sL) * K * w^2/(s^2 + 2dw s + w^2)
    // poles = -dw +- w sqrt(d^2 - 1)
    // https://www.tutorialspoint.com/control_systems/control_systems_response_second_order.htm
    pub d: f64,
    pub w: f64,
    pub K: f64,
    pub L: f64,
    pub d_lower: f64,
    pub d_upper: f64,
    pub w_lower: f64,
    pub w_upper: f64,
    pub K_lower: f64,
    pub K_upper: f64,
    pub L_lower: f64,
    pub L_upper: f64,
}

impl TransferFunction for ComplexSecondOrderSystem {
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
        let t = t - self.L;

        if t < 0.0 {
            return 0.0;
        }

        self.K * if d == 0.0 {
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

        self.K.abs() * wp.powi(2) / (((wp.powi(2) - w.powi(2)).powi(2) + (2f64*d*wp*w).powi(2)).sqrt())
    }

    fn bode_phase(&self, w: f64) -> f64 {
        use std::f64::consts::PI;

        let (d, wp) = (self.d, self.w);

        let ph = -self.L * w - (2.0*d*wp*w / (wp.powi(2) - w.powi(2))).atan();
        if ph > 0.0 {
            ph - PI
        } else {
            ph
        }
    }

    fn adjust_pole_zero(&mut self, re: f64, im: f64) {
        if re >= 0.0 {
            return
        }

        let mut d_new = self.d;
        let w_new;

        if self.d < 1.0 { // two complex poles
            let (re2, im2) = (re.powi(2), im.powi(2));
            d_new = (re2/(re2+im2)).sqrt();
            w_new = -re/d_new;
        } else if self.d == 1.0 { // real double pole
            w_new = -re;
        } else { // two real poles
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

