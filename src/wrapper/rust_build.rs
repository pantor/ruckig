fn main() {
    cxx_build::bridge("src/wrapper/rust.rs")
        .include("include")
        .include("src/wrapper")
        .include("third_party")
        .file("src/ruckig/brake.cpp")
        .file("src/ruckig/cloud_client.cpp")
        .file("src/ruckig/position_first_step1.cpp")
        .file("src/ruckig/position_first_step2.cpp")
        .file("src/ruckig/position_second_step1.cpp")
        .file("src/ruckig/position_second_step2.cpp")
        .file("src/ruckig/position_third_step1.cpp")
        .file("src/ruckig/position_third_step2.cpp")
        .file("src/ruckig/velocity_second_step1.cpp")
        .file("src/ruckig/velocity_second_step2.cpp")
        .file("src/ruckig/velocity_third_step1.cpp")
        .file("src/ruckig/velocity_third_step2.cpp")
        .flag("-DWITH_CLOUD_CLIENT=ON")
        .std("c++20")
        .compile("ruckig");

    println!("cargo:rerun-if-changed=include/ruckig");
    println!("cargo:rerun-if-changed=src/ruckig");
    println!("cargo:rerun-if-changed=src/wrapper");
}
