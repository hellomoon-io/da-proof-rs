use std::fs;
use std::io::Write;

use da_proof::gf16::GF16LUT;

fn main() {
    println!("Generating GF16 LUT...");
    let gf16_lut = GF16LUT::generate();
    println!("Writing GF16 LUT to gf16_lut.rs...");
    let mut w = fs::File::create("gf16_lut.rs").unwrap();
    write!(&mut w, "use super::gf16::{{GF16, GF16LUT}};\n").unwrap();
    write!(&mut w, "pub const LUT: GF16LUT = {:#?};", gf16_lut).unwrap();
    println!("Done.")
}
