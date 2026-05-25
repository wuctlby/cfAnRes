'''
Script to project the MC distributions and apply the pt weights from the AnRes.root of Dtask
python3 proj_thn.py config_flow.yml --cutsetConfig config_cutset.yml [-c --correlated]
If the last argument is not provided, the script will project the combined cutsets.
'''
import ROOT
import yaml
import argparse
import sys
import os
import numpy as np
from ROOT import TFile, TObject, TH1F
from alive_progress import alive_bar
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../utils")
from utils_proj import get_pt_weights, proj_multitrial, get_pt_preprocessed_sparses, proj_mc_reco, proj_mc_gen, proj_data
from utils import get_centrality_bins, make_dir_root_file, logger

ROOT.TH1.AddDirectory(False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument('--cutsetConfig', "-cc", metavar='text', type=str, nargs='?',
                        const=None, default='cutsetConfig.yaml',
                        help='Optional cutset configuration file (default: cutsetConfig.yaml)')
    parser.add_argument("--multitrial_folder", "-multfolder", metavar="text",
                        default="", help="Produce projection files for multitrial systematics")
    parser.add_argument("--multitrial_workers", "-multworkers", metavar="int",
                        type=int, default=1, help="Number of workers for multitrial projections")
    parser.add_argument("--outputDir", "-o", metavar="text",
                        default="", help="output directory, used only for directly running the script")
    args = parser.parse_args()

    # —————— Load config and operations ————————————————————
    with open(args.config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)
    operations = config["operations"]

    # —————— Check if running multitrial projections ————————————————————
    if args.multitrial_folder != "":
        logger(f"Running multitrial projections with config: {args.config}", level='INFO')
        proj_multitrial(config, args.multitrial_folder, args.multitrial_workers, operations)
        sys.exit(0)

    # —————— Get cutset info and prepare output file ————————————————————
    with open(args.cutsetConfig, 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
        iCut = f"{int(cutSetCfg['icutset']):02d}"

    outDir = os.path.join(os.path.dirname(os.path.dirname(args.cutsetConfig)), 'projs') if args.outputDir == "" else args.outputDir
    outfilePath = os.path.join(outDir, f"proj_{iCut}.root")
    os.makedirs(outDir, exist_ok=True)
    if operations.get("proj_data") or operations.get("proj_mc"):
        if os.path.exists(outfilePath):
            logger(f"Found previous projection file {outfilePath}, will update it", level='INFO')
            outfile = TFile.Open(outfilePath, 'UPDATE')
        else:
            logger(f"No previous projection file found, will create a new one at {outfilePath}", level='INFO')
            outfile = TFile(outfilePath, 'RECREATE')
    else:
        sys.exit(0)
    write_opt_data = TObject.kOverwrite if operations.get("proj_data") else 0
    write_opt_mc = TObject.kOverwrite if operations.get("proj_mc") else 0

    # —————— Get pt weights and cent di for MC projections ————————————————————
    if operations.get("proj_mc"):
        sPtWeightsD, sPtWeightsB, Bspeciesweights = get_pt_weights(config["projections"]) if config['projections'].get('PtWeightsFile') else (None, None, None)
        if config['projections'].get('CentDiffBinsMCYieldsStep'):
            proj_mc_cent_diff = True
            _, (centLowLim, centMaxLim) = get_centrality_bins(config['centrality'])

    with alive_bar(len(cutSetCfg['Pt']['min']), title='Processing pT bins') as bar:
        for i_pt, (pt_min, pt_max, bkg_min, bkg_max, fd_min, fd_max) in enumerate(zip(cutSetCfg['Pt']['min'], cutSetCfg['Pt']['max'],
                                                                                   cutSetCfg['ScoreBkg']['min'], cutSetCfg['ScoreBkg']['max'],
                                                                                   cutSetCfg['ScoreFD']['min'], cutSetCfg['ScoreFD']['max'])):

            # —————— Prepare output file and get sparses for the current pt bin ————————————————————
            logger(f'Projecting distributions for {pt_min:.1f} < pT < {pt_max:.1f} GeV/c')
            pt_label = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
            make_dir_root_file(pt_label, outfile)
            sparse_data, sparses_reco, sparses_gen, axes = get_pt_preprocessed_sparses(config, pt_label)
            outfile.cd(pt_label)

            # —————— Project data with the specified processes ————————————————————
            if operations.get("proj_data"):
                for process in config["projections"].get("process", "proj_data"): #TODO: add the possibility to specify different processes at the same time
                    proj_data(i_pt, process, sparse_data, axes, config["projections"]["proj_data"], write_opt_data)
                logger("Projected data!")

            # —————— Project MC ————————————————————————————————————————
            if operations.get("proj_mc"):
                for key, i_sparse in sparses_reco.items():
                    i_sparse.GetAxis(axes[key]['ScoreBkg']).SetRangeUser(bkg_min, bkg_max)
                    i_sparse.GetAxis(axes[key]['ScoreFD']).SetRangeUser(fd_min, fd_max)

                proj_mc_reco(sparses_reco, sPtWeightsD, sPtWeightsB, Bspeciesweights, write_opt_mc, axes, pt_min, pt_max, save_centrality=proj_mc_cent_diff)
                logger("Projected mc reco!")
                proj_mc_gen(sparses_gen, sPtWeightsD, sPtWeightsB, Bspeciesweights, write_opt_mc, axes, pt_min, pt_max, save_centrality=proj_mc_cent_diff)
                logger("Projected mc gen!\n\n")


                if proj_mc_cent_diff:
                    centStep = config['projections']['CentDiffBinsMCYieldsStep']
                    centBins = list(np.arange(centLowLim, centMaxLim + centStep, centStep))
                    hCentDiffYieldsRecoPrompt = TH1F("hCentDiffYieldsRecoPrompt", ";Centrality (%);Yield", len(centBins)-1, np.asarray(centBins, 'd'))
                    hCentDiffYieldsRecoFD = TH1F("hCentDiffYieldsRecoFD", ";Centrality (%);Yield", len(centBins)-1, np.asarray(centBins, 'd'))
                    hCentDiffYieldsGenPrompt = TH1F("hCentDiffYieldsGenPrompt", ";Centrality (%);Yield", len(centBins)-1, np.asarray(centBins, 'd'))
                    hCentDiffYieldsGenFD = TH1F("hCentDiffYieldsGenFD", ";Centrality (%);Yield", len(centBins)-1, np.asarray(centBins, 'd'))
                    for i_cent_bin, (cent_min, cent_max) in enumerate(zip(centBins[:-1], centBins[1:])):
                        cent_label = f'cent_{int(cent_min)}_{int(cent_max)}'
                        make_dir_root_file(f"{pt_label}/{cent_label}", outfile)
                        outfile.cd(f"{pt_label}/{cent_label}")
                        for key, i_sparse in sparses_reco.items():
                            i_sparse.GetAxis(axes[key]['Cent']).SetRangeUser(cent_min, cent_max)
                        for key, i_sparse in sparses_gen.items():
                            i_sparse.GetAxis(axes[key]['Cent']).SetRangeUser(cent_min, cent_max)

                        hPtPrompt, hPtFD = proj_mc_reco(sparses_reco, sPtWeightsD, sPtWeightsB, Bspeciesweights, write_opt_mc, axes, pt_min, pt_max, save_centrality=True)
                        hCentDiffYieldsRecoPrompt.SetBinContent(i_cent_bin+1, hPtPrompt.Integral())
                        hCentDiffYieldsRecoFD.SetBinContent(i_cent_bin+1, hPtFD.Integral())
                        hGenPtPrompt, hGenPtFD = proj_mc_gen(sparses_gen, sPtWeightsD, sPtWeightsB, Bspeciesweights, write_opt_mc, axes, pt_min, pt_max, save_centrality=True)
                        hCentDiffYieldsGenPrompt.SetBinContent(i_cent_bin+1, hGenPtPrompt.Integral())
                        hCentDiffYieldsGenFD.SetBinContent(i_cent_bin+1, hGenPtFD.Integral())
                        logger(f"Projected mc reco and gen for cent {cent_min}-{cent_max}!", "INFO")

                    logger("\n\n")
                    outfile.cd(pt_label)
                    hCentDiffYieldsRecoPrompt.Sumw2()
                    hCentDiffYieldsRecoPrompt.Write("hCentDiffYieldsRecoPrompt", write_opt_mc)
                    hCentDiffYieldsRecoFD.Sumw2()
                    hCentDiffYieldsRecoFD.Write("hCentDiffYieldsRecoFD", write_opt_mc)
                    hCentDiffYieldsGenPrompt.Sumw2()
                    hCentDiffYieldsGenPrompt.Write("hCentDiffYieldsGenPrompt", write_opt_mc)
                    hCentDiffYieldsGenFD.Sumw2()
                    hCentDiffYieldsGenFD.Write("hCentDiffYieldsGenFD", write_opt_mc)

                    # Restore full cent range
                    for key, i_sparse in sparses_reco.items():
                        i_sparse.GetAxis(axes[key]['Cent']).SetRangeUser(centLowLim, centMaxLim)
                    for key, i_sparse in sparses_gen.items():
                        i_sparse.GetAxis(axes[key]['Cent']).SetRangeUser(centLowLim, centMaxLim)

            bar()

    outfile.Close()
