#![warn(clippy::all, rust_2018_idioms)]
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use basic_print::basic_print; // eframe already has this kind of logging so I should probably figure out if that it expose some item for that.

// Application startin points for the two targets
#[cfg(not(target_arch = "wasm32"))]
fn main() {
    basic_print("Starting up...");

    // Log to stdout (if you run with `RUST_LOG=debug`).
    tracing_subscriber::fmt::init();

    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "eframe template",
        native_options,
        // Box::new(|cc| Box::new(control_web_apps::TemplateApp::new(cc))),
        Box::new(|cc| Box::new(control_web_apps::ControlApp::new(cc))),
    );
}

#[cfg(target_arch = "wasm32")]
fn main() {
    basic_print("Starting up...");

    // Make sure panics are logged using `console.error`.
    console_error_panic_hook::set_once();

    // Redirect tracing to console.log and friends:
    tracing_wasm::set_as_global_default();

    let web_options = eframe::WebOptions::default();
    eframe::start_web(
        "the_canvas_id", // hardcode it
        web_options,
        // Box::new(|cc| Box::new(control_web_apps::TemplateApp::new(cc))),
        Box::new(|cc| Box::new(control_web_apps::ControlApp::new(cc))),
    )
    .expect("failed to start eframe");
}
