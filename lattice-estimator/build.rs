use std::path::Path;
use std::process::Command;

fn main() {
    // TODO: propagate the PYO3_PYTHON environment variable to cargo build, which is currently not really possible in Rust
    // Option 1: write to a pyo3 config file, whose location is hard-coded in config.toml
    // Option 2: use cargo-cmake or equivalent to set the environment variable before running cargo build
    let output = Command::new("sage")
        .arg("-sh")
        .arg("-c")
        .arg("echo $SAGE_ROOT")
        .output()
        .expect("failed to retrieve SAGE_ROOT");
    let sage_python =
        Path::new(&String::from_utf8(output.stdout).unwrap().trim()).join("venv/bin/python3");

    Command::new("env")
        .arg(format!("PYO3_PYTHON={}", sage_python.to_str().unwrap()))
        .status()
        .expect("failed to export PYO3_PYTHON");

    //println!("cargo:warning=sage_python: {:?}", sage_python);
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!(
        "cargo:rustc-env=PYO3_PYTHON={}",
        sage_python.to_str().unwrap()
    );
    println!("cargo:rerun-if-env-changed=PYO3_PYTHON");
}
