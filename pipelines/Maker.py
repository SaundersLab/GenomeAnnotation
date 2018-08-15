import bioluigi.tasks
from bioluigi.decorators import requires, inherits
import luigi
import os
from luigi import LocalTarget

from . import RnaPreprocess

luigi.auto_namespace(scope=__name__)

maker_opts_1 = {'organism_type': 'eukaryotic',
                'model_org': 'fungi',
                'repeat_protein': '/tgac/software/testing/maker/2.31.8/x86_64/data/te_proteins.fasta',
                'est2genome': 1,
                'protein2genome': 1,
                'TMP': '/tmp'
                }

maker_opts_2 = {'organism_type': 'eukaryotic',
                'model_org': 'fungi',
                'repeat_protein': '/tgac/software/testing/maker/2.31.8/x86_64/data/te_proteins.fasta',
                'est2genome': 0,
                'protein2genome': 0,
                'TMP': '/tmp',
                'est_pass': 1,
                'protein_pass': 1,
                'rm_pass': 1,
                'keep_preds': 0}


class UniProt(luigi.ExternalTask):
    '''Protein evidence for MAKER'''
    def output(self):
        return LocalTarget("/usr/users/JIC_a1/buntingd/GenomeAnnotation/uniprot_PST78.fasta")


@requires(RnaPreprocess.ReferenceFasta)
class RepeatMasker(bioluigi.tasks.repeatmodeler.RepeatMasker):
    '''Repeat soft masked genome'''
    pass


# ----------- Initial MAKER run ------------- #
'''Does an intial run of maker that uses the protein and transcript evidence
   to predict genes.
   Calculates the AED distribution of this run and packages everything up for Jbrowse
'''


@inherits(RnaPreprocess.ReferenceFasta, RnaPreprocess.Mikado, UniProt, bioluigi.tasks.maker.MAKER)
class MAKER1(luigi.WrapperTask):
    maker_opts, maker_prefix = None, None

    def requires(self):
        return self.clone(requires(genome=RnaPreprocess.ReferenceFasta,
                                   est_gff=RnaPreprocess.Mikado,
                                   protein=UniProt)(bioluigi.tasks.maker.MAKER),
                          maker_opts=maker_opts_1,
                          maker_prefix="round1")

    def output(self):
        return self.requires().output()


@requires(MAKER1)
class GFFMerge1(bioluigi.tasks.maker.GFFMerge):
    pass


@requires(gff=GFFMerge1, genome=RnaPreprocess.AddGenome)
class AddGFF1(bioluigi.tasks.jbrowse.AddGFF):
    pass


@requires(GFFMerge1)
class AEDDist1(bioluigi.tasks.maker.AEDDist):
    pass


@requires(MAKER1)
class Maker2Jbrowse1(bioluigi.tasks.maker.Maker2Jbrowse):
    pass


# ----------- SNAP Training ------------- #
'''Uses the initial maker run to train SNAP ab intito gene caller
   and produce some stats.
'''


@requires(MAKER1)
class Maker2ZFF1(bioluigi.tasks.maker.Maker2ZFF):
    pass


@requires(Maker2ZFF1)
class TrainSNAP1(bioluigi.tasks.maker.TrainSNAP):
    pass


@requires(Maker2ZFF1)
class FathomStats1(bioluigi.tasks.maker.FathomStats):
    pass


# ----------- Genemark Training ------------- #
'''Genemark-ES self training'''


@requires(RepeatMasker)
class GenemarkESTrain(bioluigi.tasks.genemark.GenemarkESTrain):
    pass


# ----------- Second MAKER ------------- #
'''Do a second maker run using both the SNAP and Genemark predictors
   Calculates the AED distribution of this run and packages everything up for Jbrowse
'''


@inherits(TrainSNAP1, GenemarkESTrain)
class MAKER2(luigi.WrapperTask):
    maker_opts, maker_prefix = None, None

    def requires(self):
        return self.clone(requires(genome=RnaPreprocess.ReferenceFasta,
                                   maker_gff=GFFMerge1,
                                   snaphmm=TrainSNAP1,
                                   gmhmm=GenemarkESTrain)(bioluigi.tasks.maker.MAKER),
                          maker_opts=maker_opts_2,
                          maker_prefix="round2")

    def output(self):
        return self.requires().output()


@requires(MAKER2)
class GFFMerge2(bioluigi.tasks.maker.GFFMerge):
    pass


@requires(gff=GFFMerge2, genome=RnaPreprocess.AddGenome)
class AddGFF2(bioluigi.tasks.jbrowse.AddGFF):
    pass


@requires(GFFMerge2)
class AEDDist2(bioluigi.tasks.maker.AEDDist):
    pass


@requires(MAKER2)
class Maker2Jbrowse2(bioluigi.tasks.maker.Maker2Jbrowse):
    pass


# ----------- Second SNAP Training -------------##
'''Re-train SNAP using the output of the second maker run'''


@requires(MAKER2)
class Maker2ZFF2(bioluigi.tasks.maker.Maker2ZFF):
    pass


@requires(Maker2ZFF2)
class TrainSNAP2(bioluigi.tasks.maker.TrainSNAP):
    pass


@requires(Maker2ZFF2)
class FathomStats2(bioluigi.tasks.maker.FathomStats):
    pass


# ----------- Third MAKER -------------##
'''Do a third maker run using both the re-trained SNAP and Genemark predictors
   Calculates the AED distribution of this run and packages everything up for Jbrowse
'''


@inherits(MAKER2, TrainSNAP2)
class MAKER3(luigi.WrapperTask):
    maker_opts, maker_prefix = None, None

    def requires(self):
        return self.clone(requires(genome=RnaPreprocess.ReferenceFasta,
                                   maker_gff=GFFMerge2,
                                   snaphmm=TrainSNAP2,
                                   gmhmm=GenemarkESTrain)(bioluigi.tasks.maker.MAKER),
                          maker_opts=maker_opts_2,
                          maker_prefix="round3")

    def output(self):
        return self.requires().output()


@requires(MAKER3)
class GFFMerge3(bioluigi.tasks.maker.GFFMerge):
    pass


@requires(gff=GFFMerge3, genome=RnaPreprocess.AddGenome)
class AddGFF3(bioluigi.tasks.jbrowse.AddGFF):
    pass


@requires(MAKER3)
class Maker2Jbrowse3(bioluigi.tasks.maker.Maker2Jbrowse):
    pass


@requires(GFFMerge3)
class AEDDist3(bioluigi.tasks.maker.AEDDist):
    pass

# ----------- Pipeline Wrapper -------------##


@requires(MAKER3)
class MakerGFF3(bioluigi.tasks.maker.GFFMerge):
    '''Final output gff, include only final models'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.genes_only = True

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'maker', 'maker.gff'))


@requires(MAKER3,
          AddGFF1, AddGFF2, AddGFF3,
          FathomStats1, FathomStats2,
          AEDDist1, AEDDist2, AEDDist3, MakerGFF3)
class MAKERPipelineWrapper(luigi.WrapperTask):
    pass
