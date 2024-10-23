fn main ()
{
    cxx_build::bridge("src/lib.rs")
        .file("src/SerialDeserial.cc")
        .include("/usr/local/include/openfhe")
        .include("/usr/local/include/openfhe/third-party/include")
        .include("/usr/local/include/openfhe/core")
        .include("/usr/local/include/openfhe/pke")
        .include("/usr/local/include/openfhe/binfhe")
        .include("./openfhe-development/install/include") // GitHub Actions
        .include("./openfhe-development/install/include/openfhe") // GitHub Actions
        .include("./openfhe-development/install/include/openfhe/third-party/include") // GitHub Actions
        .include("./openfhe-development/install/include/openfhe/core") // GitHub Actions
        .include("./openfhe-development/install/include/openfhe/pke") // GitHub Actions
        .include("./openfhe-development/install/include/openfhe/binfhe") // GitHub Actions
        .flag_if_supported("-std=c++17")
        .flag_if_supported("-Wall")
        .flag_if_supported("-Werror")
        .flag_if_supported("-O3")
        .flag_if_supported("-fopenmp") // [-Wunknown-pragmas]
        .flag_if_supported("-Wno-parentheses") // [-Wparentheses]
        .flag_if_supported("-Wno-unused-parameter") // [-Wunused-parameter]
        .flag_if_supported("-Wno-missing-field-initializers") // [-Wmissing-field-initializers]
        .flag_if_supported("-Wno-unused-function") // [-Wunused-function]
        .compile("openfhe");

    // linking openFHE
    println!("cargo::rustc-link-arg=-L/usr/local/lib");
    println!("cargo::rustc-link-arg=-lOPENFHEpke");
    println!("cargo::rustc-link-arg=-lOPENFHEbinfhe");
    println!("cargo::rustc-link-arg=-lOPENFHEcore");
    // linking OpenMP
    println!("cargo::rustc-link-arg=-fopenmp");
    // necessary to avoid LD_LIBRARY_PATH
    println!("cargo::rustc-link-arg=-Wl,-rpath=/usr/local/lib");
}