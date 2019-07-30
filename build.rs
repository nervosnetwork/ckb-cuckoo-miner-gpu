extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/worker/cuckoo.c")
        .include("src/worker")
        .flag("-O3")
        .flag("-lcrypto")
        // .flag("-mavx")
        // .flag("-mavx2")
        // .flag("-march=native")
        // .flag("-mavx512f")
        // .flag("-mavx512cd")
        .static_flag(true)
        .compile("libmean.a");
}