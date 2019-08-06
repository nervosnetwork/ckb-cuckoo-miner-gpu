extern crate cc;

fn main() {
    cc::Build::new()
        .file("src/worker/cuckoo.cu")
        .include("src/worker")
        .flag("-O3")
        .flag("-lcrypto")
        .cuda(true)
        .compile("libmean.a");

    // Add link directory
    // - This path depends on where you install CUDA (i.e. depends on your Linux distribution)
    // - This should be set by `$LIBRARY_PATH`
    println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    println!("cargo:rustc-link-lib=cudart");
}
