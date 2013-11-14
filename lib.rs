#[feature(macro_rules)];

#[link(name="bzip2",
    vers="0.1",
    url="https://github.com/MarkJr94/bzip2rs",
    package_id="bzip2rs")];

#[comment = "Rust port of BZIP2 algorithm"];
#[license = "MIT"];
#[crate_type = "lib"];

// #[deny(non_camel_case_types)];
#[deny(non_uppercase_statics)];
#[deny(unnecessary_qualification)];
#[warn(missing_doc)];

pub use self::bzip2reader::Bzip2Reader;

mod macros;
mod consts;
mod checksum;
mod bzip2reader;

