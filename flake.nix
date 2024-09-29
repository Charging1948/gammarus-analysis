{
  description = "Development environment for data_analysis_revamp.R";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
  };
  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
      r_pkgs = with pkgs.rPackages; [
          ggplot2
          dplyr
          tidyverse
          readxl
          janitor
          ggthemes
          lubridate
          GGally
          survminer
          car
          gridExtra
          gridtext
          multcomp
          ggsurvfit
          gtsummary
          condSURV
          tidycmprsk
          flexsurv
          drc
      ];
      rstudio-wrapped = pkgs.rstudioWrapper.override {
        packages = r_pkgs;
      };
    in
    {
      devShell.x86_64-linux = pkgs.mkShell {
        buildInputs = with pkgs.rPackages; [
          rstudio-wrapped
          
          pkgs.R
        ] ++ r_pkgs;

        shellHook = ''
          echo "Entering the R development environment..."
        '';
      };
    };
}
