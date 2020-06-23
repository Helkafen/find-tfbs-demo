{ pkgs }:
with pkgs;
rustPlatform.buildRustPackage rec {
    pname = "find-tfbs";
    version = "1.0.1";

    nativeBuildInputs = [ pkgconfig clang ];
    buildInputs = [
      llvmPackages.libclang
      zlib
    ];

    # Verify the packages that find-tfbs depend on
    cargoSha256 = "1k7hgmp5598czc492wb3nfgqqfbqw8hc43h4gpqwm83fib3rf877";
    verifyCargoDeps = true;

    src = builtins.fetchGit {
      url = https://github.com/Helkafen/find-tfbs;
      rev = "613ddfd1be4616c57256bab17a4a8e7d263627cf";
    };

    doCheck = false; # Don't run the tests here. They depend on HOCOMOCO files that are not distributed with the find-tfbs source code
    
    RUST_BACKTRACE = 1; # In case of a runtime error, show a full stack trace
    LIBCLANG_PATH = llvmPackages.libclang + "/lib";
  }
