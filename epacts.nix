{ pkgs }:
with pkgs;
stdenv.mkDerivation {
  name = "epacts";

  nativeBuildInputs = [ autoreconfHook ];
  buildInputs = [
    gcc
    gmp
    zlib
    bzip2
    lzma
    autoreconfHook
    automake
    perl
    ghostscript
    groff
    gnuplot
    R
  ];

  src = builtins.fetchGit {
    url = https://github.com/statgen/EPACTS;
    rev = "edc39656dec3369a168315bfee128d617dde1135";
  };

  prePatch = ''
        sed -i 's/1.14/1.15/g' configure;
        sed -i 's/--foreign/--foreign --add-missing --force-missing/g' Makefile;
        sed -i 's/--foreign/--foreign --add-missing --force-missing/g' Makefile.in
  '';

}

