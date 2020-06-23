# find-tfbs use case

This project shows how to use find-tfbs to find associations between blood traits and variation in TFBS counts. The search is performed over the open chromatin regions of several hematopoietic cell types.

## How to install the required packages

### Option 1: Install the Nix package manager, then run:
```console
$ nix-shell
```

This will download and install all requirements automatically and enter a new shell where these packages are available. This method is perfectly reproducible.

## Option 2: Install the packages manually. The requirements for EPACTS and find-tfbs are described in epacts.nix and find-tfbs.nix
```console
$ git clone https://github.com/statgen/EPACTS.git
$ cd EPACTS/; 
$ sed -i 's/1.14/1.15/g' configure;
sed -i 's/--foreign/--foreign --add-missing --force-missing/g' Makefile;
sed -i 's/--foreign/--foreign --add-missing --force-missing/g' Makefile.in;
$ ./configure --prefix `pwd`/epacts
$ make -j 4
$ make install -j 4
```

```console
$ git clone git@github.com:Helkafen/find-tfbs.git
$ cd find-tfbs
$ git ckeckout fd728386872b28f64308873158a4c32f41ef7e57
$ ./build_static.sh
$ cp target/x86_64-unknown-linux-musl/release/find-tfbs ..
$ cd ..
```

## How to run the pipeline

### On a SLURM cluster:
```console
snakemake --cluster sbatch -A <cluster account id> -t {resources.runtime} --mem-per-cpu {resources.mem} -J TFBS --jobs 4 all
```

### In the current shell session:
```console
snakemake --jobs 4 all
```