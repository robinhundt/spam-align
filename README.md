# SpaM-Align

![Rust](https://github.com/robinhundt/spam-align/workflows/Rust/badge.svg)

In order to compile this program, the latest stable Rust toolchain is needed. Installing this can be conveniently done via the toolchain manager `rustup` [(Link)](https://rustup.rs/).  

After having installed the latest stable toolchain, the program can be built and run by executing
```
cargo run --release -- --help
```

Note that program specific parameters must be passed after the `--`.  
Alternatively the binary (called `spam-align`) can be installed by issuing:
```
cargo install --path .
```

## Tests

The contained tests can be run by executing:
```
cargo test --release
```
Note that, while `--release` is not required, the tests take considerably longer if executed in debug mode.


## Benchmarks

There are some benchmarks located under `benches` which can be run with:
```
cargo bench
``` 

## Logging

When setting the Env variable `RUST_LOG=info` the program will print some timing information.


## Notes
The algorithm for checking and updating transitivity frontiers is a reimplementation and improvement ofgst
 [GABIOS-LIB](gobics.de/burkhard/papers/jobim.pdf).  
