let
  pkgs = import (fetchTarball https://github.com/NixOS/nixpkgs/archive/17.03.tar.gz) {};
  pkgs3 = import (fetchTarball https://github.com/NixOS/nixpkgs/archive/20.03.tar.gz) {};
  tabix = import ./tabix.nix { pkgs=pkgs3; };
  find-tfbs = import ./find-tfbs.nix { pkgs=pkgs3; };
  epacts = import ./epacts.nix { pkgs=pkgs; };
in with pkgs; {
  simpleEnv = stdenv.mkDerivation {
    name = "tfbs-env";
    version = "1";
    buildInputs = [
      pkgs.R
      pkgs.rPackages.ggplot2
      pkgs.rPackages.tidyverse
      pkgs.gnuplot
      epacts
      tabix
      find-tfbs
      pkgs3.vcftools
      pkgs3.bcftools
      pkgs3.snakemake
      pkgs3.samtools
      pkgs.python3
      pkgs3.bedtools
    ];

  };

  LANG="C";
}
