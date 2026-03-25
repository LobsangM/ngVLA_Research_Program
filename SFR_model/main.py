from rebin_beam_sampling import SFRRadioModel

if __name__ == "__main__":

    model = SFRRadioModel("config.yaml")

    model.luminosity_to_flux()

    model.plot_sfr()
    
    model.plot_flux()

    model.print_statistics()

    model.save_flux_to_fits()