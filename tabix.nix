{ pkgs }:
with pkgs;
stdenv.mkDerivation {
    name = "tabix";
    nativeBuildInputs = [ autoreconfHook ];
    buildInputs = [
      gcc
      zlib
      bzip2
      lzma
      autoreconfHook
    ];

    src = builtins.fetchGit {
      url = https://github.com/samtools/htslib;
      rev = "2b6bac81b45af266b8cf3e851df08eb4f53c7157";
    };
    
  }
