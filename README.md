# Spam Align

In order to compile this program, the latest stable Rust toolchain is needed. Installing this can be conveniently done via the toolchain manager `rustup` [(Link)](https://rustup.rs/).  

After having installed the latest stable toolchain, the program can be built and run by executing
```
cargo run --release --bin stats -- --help
```

Note that program specific parameters must be passed after the `--`.

In order to use the program the `bb3_release` folder of the BAliBASE version 3 dataset http://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz should be placed in a `data` directory within the root of this project. 

## Notes
The algorithm for checking and updating consistency bound is based on [GABIOS-LIB](gobics.de/burkhard/papers/jobim.pdf).  


## Benchmakrk data todo
OxBase