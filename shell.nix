{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    python3
    python311Packages.networkx
		python311Packages.matplotlib
    fish
  ];

  shellHook = ''
    exec fish
  '';
}
