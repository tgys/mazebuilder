# shell.nix
{ pkgs ? import <nixpkgs> {} }:
let
  my-python = pkgs.python2;
#  glad = pkgs.callPackage /home/tei/projects/nixpkgs/pkgs/development/python-modules/glad/default.nix {};
  python-with-my-packages = my-python.withPackages (p: with p; [
    urllib3
    # other python packages you want
    #
  ]);
in
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.gcc
    pkgs.boost-build
    pkgs.cmake
  ];
  buildInputs = [
    python-with-my-packages
    pkgs.libGLU
    pkgs.minizip
    pkgs.assimp
    pkgs.freetype
    pkgs.xorg.libXext
    pkgs.doxygen
    pkgs.libglvnd
    pkgs.boost
    pkgs.libGL
    pkgs.xorg.libX11.dev
    pkgs.mesa.dev
    pkgs.xorg.libXrandr
    pkgs.xorg.libXinerama
    pkgs.xorg.libXcursor
    pkgs.xorg.xinput
    pkgs.xorg.libXi
    pkgs.xorg.libXfixes
    pkgs.freeglut
    pkgs.glfw
    pkgs.stb
  ];
  shellHook = ''
    PYTHONPATH=${python-with-my-packages}/${python-with-my-packages.sitePackages}
    # maybe set more env-vars
  '';
}
