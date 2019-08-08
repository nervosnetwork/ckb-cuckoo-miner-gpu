fn main() {
    // Add link directory
    // - This path depends on where you install CUDA (i.e. depends on your Linux distribution)
    // - This should be set by `$LIBRARY_PATH`
    println!("cargo:rustc-link-search=native=C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v10.1\\lib\\x64");
    println!("cargo:rustc-link-lib=cudart");
    println!("cargo:rustc-link-search=native=miner\\src\\worker");
}
