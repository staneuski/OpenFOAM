{
  description = "Local OpenFOAM environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-darwin";
      pkgs = import nixpkgs { inherit system; };
    in {
      devShells.${system}.default = (import nixpkgs { inherit system; }).mkShell {
        packages = [
          pkgs.boost
          pkgs.cgal
          pkgs.cmake
          pkgs.metis
          pkgs.mpi
          pkgs.scotch
          pkgs.zlib
          pkgs.flex
        ];

        shellHook = ''
          [ -f etc/bashrc ] &&
            source etc/bashrc
        '';
      };
    };
}
