#![warn(clippy::all, rust_2018_idioms)]

mod transfer_functions;

#[allow(unused_imports)]
use basic_print::basic_print; // basic print for print-debugging

use frequency_response_app::FreqResp;
use pole_position_app::PolePos;

pub struct ControlApp {
    cur_app_idx: Option<usize>,
    apps: Vec<Box<dyn CentralApp>>,
}

trait CentralApp {
    fn get_label(&self) -> &str;
    fn draw_app(&mut self, ui: &mut egui::Ui);
}

impl eframe::App for ControlApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Set light theme so we don't have to mess with colors for plots
        ctx.set_visuals(egui::Visuals::light());

        egui::TopBottomPanel::top("app_selection_panel").show(ctx, |ui| {
            if self.top_bar(ui) {
                #[cfg(not(target_arch = "wasm32"))] // no quit on web pages!
                _frame.close();
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.centered_and_justified(|ui| {
                match self.cur_app_idx {
                    None => {
                        ui.label("Select an application in the bar above.");
                    }
                    Some(idx) => {
                        self.apps[idx].draw_app(ui);
                    }
                };
            });
        });
    }
}

impl ControlApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        // This is also where you can customized the look at feel of egui using
        // `cc.egui_ctx.set_visuals` and `cc.egui_ctx.set_fonts`.

        let apps: Vec<Box<dyn CentralApp>> = vec![
            Box::new(PolePos::new("Pole Positioning".to_string())),
            Box::new(FreqResp::new("Frequency Response".to_string())),
        ];

        if cfg!(debug_assertions) {
            ControlApp {
                cur_app_idx: Some(0),
                apps,
            }
        } else {
            ControlApp {
                cur_app_idx: None,
                apps,
            }
        }
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) -> bool {
        #[allow(unused_mut)]
        let mut quit = false;

        ui.horizontal_wrapped(|ui| {
            ui.heading("Control Apps");
            ui.separator();

            for (idx, app) in self.apps.iter().enumerate() {
                let checked = match self.cur_app_idx {
                    Some(cur) => cur == idx,
                    None => false,
                };
                if ui.selectable_label(checked, app.get_label()).clicked() {
                    if checked {
                        self.cur_app_idx = None;
                    } else {
                        self.cur_app_idx = Some(idx);
                    }
                }
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                #[cfg(not(target_arch = "wasm32"))] // no quit on web pages!
                {
                    if ui.button("Quit").clicked() {
                        quit = true;
                    }
                    ui.separator();
                }
                egui::warn_if_debug_build(ui);
            });
        });

        return quit;
    }
}

mod pole_position_app {
    #[allow(unused_imports)]
    use basic_print::basic_print; // basic print for print-debugging

    use egui::{Ui, Vec2};

    use crate::transfer_functions::*;
    use crate::CentralApp;

    use super::tf_plots;

    #[derive(PartialEq, Debug, Clone, Copy)]
    enum Systems {
        FirstOrder,
        RealSecondOrder,
        RealSecondOrderAndZero,
        ComplexSecondOrder,
    }

    #[derive(PartialEq, Debug, Clone, Copy)]
    enum Display {
        StepResponse,
        BodeDiagram,
    }

    #[derive(Debug)]
    pub struct PolePos {
        label: String,

        systype: Systems,
        display: Display,

        first: FirstOrderSystem,
        second: RealSecondOrderSystem,
        second_zero: RealSecondOrderSystemAndZero,
        complex: ComplexSecondOrderSystem,

        pole_drag_offset: Option<(f64, f64)>,
    }

    impl PolePos {
        pub fn new(label: String) -> PolePos {
            PolePos {
                label,
                systype: Systems::FirstOrder,
                display: Display::StepResponse,
                first: FirstOrderSystem { 
                    T: 1.0, T_lower: 0.1, T_upper: 500.0,
                    K: 1.0, K_lower: -2.0, K_upper: 2.0,
                    L: 0.0, L_lower: 0.0, L_upper: 5.0,
                },
                second: RealSecondOrderSystem {
                    T1: 1.0, T2: 2.0, T_lower: 0.1, T_upper: 500.0,
                    K: 1.0, K_lower: -2.0, K_upper: 2.0,
                    L: 0.0, L_lower: 0.0, L_upper: 5.0,
                },
                second_zero: RealSecondOrderSystemAndZero {
                    T1: 0.5, T2: 0.75, T_lower: 0.1, T_upper: 500.0,
                    Tz: 1.0, Tz_lower: -10.0, Tz_upper: 10.0,
                    K: 1.0, K_lower: -2.0, K_upper: 2.0,
                    L: 0.0, L_lower: 0.0, L_upper: 5.0,
                },
                complex: ComplexSecondOrderSystem {
                    d: 0.7, d_lower: 0.01, d_upper: 5.0,
                    w: 1.0, w_lower: 0.01, w_upper: 5.0,
                    K: 1.0, K_lower: -2.0, K_upper: 2.0,
                    L: 0.0, L_lower: 0.0, L_upper: 5.0,
                },
                pole_drag_offset: None,
            }
        }

        fn pole_plot(&mut self, ui: &mut Ui, width: f32, height: f32) {
            let (dragged, pointer_coordinate) = match self.systype {
                Systems::FirstOrder => tf_plots::pole_plot(&self.first, ui, width, height),
                Systems::RealSecondOrder => tf_plots::pole_plot(&self.second, ui, width, height),
                Systems::RealSecondOrderAndZero => tf_plots::pole_plot(&self.second_zero, ui, width, height),
                Systems::ComplexSecondOrder => tf_plots::pole_plot(&self.complex, ui, width, height),
            };

            // Handle dragging
            if dragged {
                if let Some((re,im)) = pointer_coordinate { // This should never fail
                    match self.systype {
                        Systems::FirstOrder => self.first.adjust_pole_zero(re, im),
                        Systems::RealSecondOrder => self.second.adjust_pole_zero(re, im),
                        Systems::RealSecondOrderAndZero => self.second_zero.adjust_pole_zero(re, im),
                        Systems::ComplexSecondOrder => self.complex.adjust_pole_zero(re, im),
                    };
                }
            } else {
                self.pole_drag_offset = None;
            }
        }

        fn step_response_plot(&mut self, ui: &mut Ui, width: f32, height: f32) {
            let (_dragged, _pointer_coordinate) = match self.systype {
                Systems::FirstOrder => tf_plots::step_response_plot(&self.first, ui, width, height),
                Systems::RealSecondOrder => tf_plots::step_response_plot(&self.second, ui, width, height),
                Systems::RealSecondOrderAndZero => tf_plots::step_response_plot(&self.second_zero, ui, width, height),
                Systems::ComplexSecondOrder => tf_plots::step_response_plot(&self.complex, ui, width, height),
            };
        }

        fn bode_plot(&mut self, ui: &mut Ui, width: f32, height: f32) {
            let (_amp_dragged, _amp_pointer, _ph_dragged, _ph_pointer) = match self.systype {
                Systems::FirstOrder => tf_plots::bode_plot(&self.first, ui, width, height),
                Systems::RealSecondOrder => tf_plots::bode_plot(&self.second, ui, width, height),
                Systems::RealSecondOrderAndZero => tf_plots::bode_plot(&self.second_zero, ui, width, height),
                Systems::ComplexSecondOrder => tf_plots::bode_plot(&self.complex, ui, width, height),
            };
        }

        fn order_selection(&mut self, ui: &mut Ui) {
            egui::ComboBox::from_label("System type")
                .selected_text(match self.systype {
                    Systems::FirstOrder => "First order",
                    Systems::RealSecondOrder => "Second order",
                    Systems::RealSecondOrderAndZero => "Second order with zero",
                    Systems::ComplexSecondOrder => "Second order complex",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.systype, Systems::FirstOrder, "First order");
                    ui.selectable_value(&mut self.systype, Systems::RealSecondOrder, "Second order");
                    ui.selectable_value(&mut self.systype, Systems::RealSecondOrderAndZero, "Second order with zero");
                    ui.selectable_value(&mut self.systype, Systems::ComplexSecondOrder, "Second order complex");
                }
            );
        }

        fn display_selection(&mut self, ui: &mut Ui) {
            ui.heading("Select Display");
            ui.horizontal(|ui| {
                ui.radio_value(&mut self.display, Display::StepResponse, "Step Response");
                ui.radio_value(&mut self.display, Display::BodeDiagram, "Bode Diagram");
            });
        }

        fn parameter_sliders(&mut self, ui: &mut Ui) {
            match self.systype {
                Systems::FirstOrder => {
                    ui.heading("G(s) = K * exp(-Ls) / (s*T + 1)");
                    ui.add(egui::Slider::new(&mut self.first.T, self.first.T_lower..=self.first.T_upper).text("T").logarithmic(true));
                    ui.add(egui::Slider::new(&mut self.first.K, self.first.K_lower..=self.first.K_upper).text("K"));
                    ui.add(egui::Slider::new(&mut self.first.L, self.first.L_lower..=self.first.L_upper).text("L"));
                }
                Systems::RealSecondOrder => {
                    ui.heading("G(s) = K * exp(-Ls) / ((s*T1 + 1) * (s*T2 + 1))");
                    ui.add(egui::Slider::new(&mut self.second.T1, self.second.T_lower..=self.second.T_upper).text("T1").logarithmic(true));
                    ui.add(egui::Slider::new(&mut self.second.T2, self.second.T_lower..=self.second.T_upper).text("T2").logarithmic(true));
                    ui.add(egui::Slider::new(&mut self.second.K, self.second.K_lower..=self.second.K_upper).text("K"));
                    ui.add(egui::Slider::new(&mut self.second.L, self.second.L_lower..=self.second.L_upper).text("L"));
                }
                Systems::RealSecondOrderAndZero => {
                    ui.heading("G(s) = K * exp(-Ls) * (s*Tz-1) / ((s*T1 + 1) * (s*T2 + 1))");
                    ui.add(egui::Slider::new(&mut self.second_zero.T1, self.second_zero.T_lower..=self.second_zero.T_upper).text("T1").logarithmic(true));
                    ui.add(egui::Slider::new(&mut self.second_zero.T2, self.second_zero.T_lower..=self.second_zero.T_upper).text("T2").logarithmic(true));
                    ui.add(egui::Slider::new(&mut self.second_zero.Tz, self.second_zero.Tz_lower..=self.second_zero.Tz_upper).text("Tz"));
                    ui.add(egui::Slider::new(&mut self.second_zero.K, self.second_zero.K_lower..=self.second_zero.K_upper).text("K"));
                    ui.add(egui::Slider::new(&mut self.second_zero.L, self.second_zero.L_lower..=self.second_zero.L_upper).text("L"));
                }
                Systems::ComplexSecondOrder => {
                    ui.heading("G(s) = K * exp(-Ls) * ω^2 / (s^2 + 2δωs + ω^2)");
                    ui.add(egui::Slider::new(&mut self.complex.d, self.complex.d_lower..=self.complex.d_upper).text("δ"));
                    ui.add(egui::Slider::new(&mut self.complex.w, self.complex.w_lower..=self.complex.w_upper).text("ω"));
                    ui.add(egui::Slider::new(&mut self.complex.K, self.complex.K_lower..=self.complex.K_upper).text("K"));
                    ui.add(egui::Slider::new(&mut self.complex.L, self.complex.L_lower..=self.complex.L_upper).text("L"));
                }
            };
        }
    }

    impl CentralApp for PolePos {
        fn draw_app(&mut self, ui: &mut Ui) {
            // ui.spacing_mut().item_spacing.x = 10.0;

            let max_width = 550.0;
            let Vec2 { x, y } = ui.available_size();
            let is_vertical = x < max_width;

            if is_vertical {
                egui::Grid::new("app_grid").num_columns(1).show(ui, |ui| {
                    ui.vertical(|ui| {
                        self.order_selection(ui);
                        ui.separator();
                        self.display_selection(ui);
                        ui.separator();
                        self.parameter_sliders(ui);
                        ui.separator();
                    });
                    ui.end_row();

                    let Vec2 { x, y } = ui.available_size();
                    let mut width = x;
                    let mut height = y / 2.0;

                    if width >= height * 1.75 {
                        width = height * 1.75
                    } else {
                        height = width / 1.75
                    }

                    self.pole_plot(ui, width, height);
                    ui.end_row();

                    match self.display {
                        Display::StepResponse => self.step_response_plot(ui, width, height),
                        Display::BodeDiagram => self.bode_plot(ui, width, height),
                    };
                });
            } else {
                let mut width = (x / 2.0).min(max_width);
                let mut height = y / 2.0;

                if width >= height * 1.75 {
                    width = height * 1.75
                } else {
                    height = width / 1.75
                }

                egui::Grid::new("app_grid").num_columns(2).show(ui, |ui| {
                    ui.vertical(|ui| {
                        self.order_selection(ui);
                        ui.add_space(20.0);
                        self.parameter_sliders(ui);
                    });

                    self.pole_plot(ui, width, height);

                    ui.end_row();
                    self.step_response_plot(ui, width, height);
                    self.bode_plot(ui, width, height);
                });
            }
        }

        fn get_label(&self) -> &str {
            &self.label
        }
    }
}




mod tf_plots {
    use egui::plot::{ Line, LineStyle, MarkerShape, Plot, PlotPoints, PlotUi, Points, };
    use egui::{ Align, Color32, InnerResponse, Layout, Ui, Vec2, };

    use std::f64::consts::PI;
    use std::ops::Range;

    use crate::transfer_functions::*;

    // Helper that give a sane default plot window. Looks can be modified with the second to last
    // argument and what is plotted is given by the last. Returns whether the plot is dragged by
    // the mouse and the plot coordinate of the mouse.
    fn plot_show(
        ui: &mut Ui,
        title: &str,
        width: f32,
        height: f32,
        x_bounds: Range<f64>,
        y_bounds: Range<f64>,
        plot_mod_fn: impl FnOnce(Plot) -> Plot,
        build_fn: impl FnOnce(&mut PlotUi),
    ) -> (bool, Option<(f64, f64)>)
    {
        let InnerResponse {
            response: _,
            inner: (dragged, pointer_coordinate),
        } = ui.allocate_ui_with_layout(
            Vec2 {
                x: width,
                y: height,
            },
            Layout::top_down(Align::LEFT),
            |ui| {
                ui.heading(title);

                let mut plot = Plot::new(title)
                    .allow_scroll(false)
                    .allow_zoom(false)
                    .allow_boxed_zoom(false)
                    .allow_drag(false)
                    .show_x(false)
                    .show_y(false)
                    .include_x(x_bounds.start)
                    .include_x(x_bounds.end)
                    .include_y(y_bounds.start)
                    .include_y(y_bounds.end)
                    .set_margin_fraction(Vec2 { x: 0.0, y: 0.0 });

                plot = plot_mod_fn(plot);

                let InnerResponse {
                    response: show_response,
                    inner: pointer_coordinate,
                } = plot.show(ui, |plot_ui| {
                    build_fn(plot_ui);
                    plot_ui.pointer_coordinate().map(|pp| (pp.x, pp.y))
                });

                (show_response.dragged(), pointer_coordinate)
            },
        );

        (dragged, pointer_coordinate)
    }

    pub fn pole_plot(
        tf: &impl TransferFunction,
        ui: &mut Ui,
        width: f32,
        height: f32,
    ) -> (bool, Option<(f64, f64)>)
    {
        // Plot params
        let marker_radius = 10.0;
        let re_bounds = -3.55..1.1;
        let im_bounds = -1.5..1.5;

        // Plot points
        let pole_data = Points::new(tf.poles());
        let zero_data = Points::new(tf.zeros());
        let unit_circle = Line::new(PlotPoints::from_parametric_callback(
            |t| (t.sin(), t.cos()),
            0.0..(2.1 * PI),
            100,
        ));

        // Plot
        plot_show(
            ui,
            "Pole Placement",
            width,
            height,
            re_bounds,
            im_bounds,
            |plot| plot.data_aspect(1.0),
            |plot_ui| {
                plot_ui.line(unit_circle.color(Color32::GRAY));
                plot_ui.points(
                    pole_data.shape(MarkerShape::Cross)
                        .color(Color32::BLACK)
                        .radius(marker_radius),
                );
                plot_ui.points(
                    zero_data.shape(MarkerShape::Circle)
                        .color(Color32::BLACK)
                        .radius(marker_radius)
                        .filled(false),
                );
            },
        )
    }

    pub fn step_response_plot(
        tf: &impl TransferFunction,
        ui: &mut Ui,
        width: f32,
        height: f32,
    ) -> (bool, Option<(f64, f64)>)
    {
        // Plot params
        let n_samples = 100;
        let t_end = 10.0;
        let pad_ratio = 0.1;

        // Calculate plot bounds
        let t_bounds = (0.0 - t_end * pad_ratio)..(t_end + t_end * pad_ratio);
        let y_bounds = (0.0 - pad_ratio)..(1.5 + pad_ratio);

        // Calc plot data
        let step = (t_bounds.end - t_bounds.start) / ((n_samples - 1) as f64);
        let mut points: Vec<[f64; 2]> = Vec::new();
        for i in 0..n_samples {
            let t = t_bounds.start + step * (i as f64);
            points.push([t, tf.step_response(t)]);
        }
        let data = Line::new(points);

        // Plot
        plot_show(
            ui,
            "Step Response",
            width,
            height,
            t_bounds,
            y_bounds,
            |plot| plot,
            |plot_ui| {
                plot_ui.line(data.color(Color32::RED).style(LineStyle::Solid));
            },
        )
    }

    pub fn bode_plot(
        tf: &impl TransferFunction,
        ui: &mut Ui,
        width: f32,
        height: f32,
    ) -> (bool, Option<(f64, f64)>, bool, Option<(f64, f64)>)
    {
        // Plot params
        let n_samples = 200;
        let w_bounds_exp = -3.0..2.0;

        // Calc plot data
        let step = (w_bounds_exp.end - w_bounds_exp.start) / ((n_samples - 1) as f64);
        let mut amp_points: Vec<[f64; 2]> = Vec::new();
        let mut phase_points: Vec<[f64; 2]> = Vec::new();
        for i in 0..n_samples {
            let we = w_bounds_exp.start + step * (i as f64);
            let w = 10f64.powf(we);
            amp_points.push([we, tf.bode_amplitude(w).log10()]);
            phase_points.push([we, tf.bode_phase(w)]);
        }
        let amp_data = Line::new(amp_points);
        let phase_data = Line::new(phase_points);

        // Plot
        let InnerResponse {
            response: _,
            inner: (amp_dragged, amp_pointer, ph_dragged, ph_pointer),
        } = ui.allocate_ui_with_layout(
            Vec2 {
                x: width,
                y: height,
            },
            Layout::top_down(Align::LEFT),
            |ui| {
                let height = (height - ui.spacing().item_spacing.y) / 2.0;
                let (amp_dragged, amp_pointer) = plot_show(
                    ui,
                    "Bode Plot - Amplitude",
                    width,
                    height,
                    w_bounds_exp.clone(),
                    -4.0..15f64.log10(),
                    |plot| plot,
                    |plot_ui| {
                        plot_ui.line(amp_data.color(Color32::RED).style(LineStyle::Solid));
                    },
                    );
                let (ph_dragged, ph_pointer) = plot_show(
                    ui,
                    "Bode Plot - Phase",
                    width,
                    height,
                    w_bounds_exp.clone(),
                    -PI / 0.95..PI / 4.0,
                    |plot| plot,
                    |plot_ui| {
                        plot_ui.line(phase_data.color(Color32::RED).style(LineStyle::Solid));
                    },
                    );
                (amp_dragged, amp_pointer, ph_dragged, ph_pointer)
            },
            );

        (amp_dragged, amp_pointer, ph_dragged, ph_pointer)
    }
}





mod frequency_response_app {
    use crate::CentralApp;

    pub struct FreqResp {
        label: String,
    }

    impl FreqResp {
        pub fn new(label: String) -> FreqResp {
            FreqResp { label }
        }
    }

    impl CentralApp for FreqResp {
        fn draw_app(&mut self, ui: &mut egui::Ui) {
            ui.label("Frequency response app currently not implemented.");
        }

        fn get_label(&self) -> &str {
            &self.label
        }
    }
}
