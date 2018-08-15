import bioluigi.tasks
from bioluigi.decorators import requires
import os
import luigi
from luigi import LocalTarget, Parameter

from . import RnaPreprocess

luigi.auto_namespace(scope=__name__)


@requires(RnaPreprocess.ReferenceFasta)
class RepeatMasker(bioluigi.tasks.repeatmodeler.RepeatMasker):
    pass


@requires(fasta=RepeatMasker, gff=RnaPreprocess.MikadoGFF3)
class CodingQuarry(bioluigi.tasks.codingquarry.CodingQuarry):
    base_dir = Parameter()

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'CodingQuarry'))


@requires(fasta=RepeatMasker, gff=RnaPreprocess.MikadoGFF3)
class CodingQuarryPM(bioluigi.tasks.codingquarry.CodingQuarryPM):
    base_dir = Parameter()

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'CodingQuarry'))


@requires(CodingQuarryPM)
class ConcatGFF(bioluigi.tasks.codingquarry.ConcatGFF):
    pass


@requires(ConcatGFF)
class CQGFF(bioluigi.tasks.codingquarry.FixGFF):
    pass


@requires(gff=CQGFF, genome=RnaPreprocess.AddGenome)
class AddGFF(bioluigi.tasks.jbrowse.AddGFF):
    pass


@requires(CQGFF, AddGFF)
class CodingQuarryWrapper(luigi.WrapperTask):
    pass
