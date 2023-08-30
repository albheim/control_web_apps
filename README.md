# Control Web Applications

A simple web app providing interactive tools to explore basic concepts in control theory for the basic course in control given at Lund University. Try it [here](https://albheim.github.io/control_web_apps/).

# Development

The app is written in Rust using [egui](https://github.com/emilk/egui/) and [eframe](https://github.com/emilk/egui/tree/master/crates/eframe). It is based on the [eframe_template](https://github.com/emilk/eframe_template/tree/master).

## Compilation - WASM

[Trunk](https://trunkrs.dev/) is used to compile the rust code to Wasm and then package the necessary HTML and JavaScript wrappers into a complete webpage. Trunk can be intalled via cargo, i.e., `cargo install --locked trunk`. Make sure the `wasm32-unknown-unknown` target for `rustc` is installed, if you are using `rustup` this can be done with `rustup target add wasm32-unknown-unknown`.

* `trunk serve` will serve the webpage on `127.0.0.1:8080` and automatically rebuild if any files are changed.
* `trunk build --release` will build a release version of the webpage into the `dist` directory.

Note, the JS wrapper is set up to cache the Wasm app which cause problem when developing. The caching can be bypassed by requesting the `index.html#dev` page.


## Compilation - Local

The [eframe](https://github.com/emilk/egui/tree/master/crates/eframe) framework also supports being compiled down to an ordinary desktop app. This can be done via normal cargo commands, e.g., `cargo run`, `cargo build --release` etc.

However, on Linux you need to first run:

`sudo apt-get install libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libspeechd-dev libxkbcommon-dev libssl-dev`

On Fedora Rawhide you need to run:

`dnf install clang clang-devel clang-tools-extra speech-dispatcher-devel libxkbcommon-devel pkg-config openssl-devel libxcb-devel fontconfig-devel`
