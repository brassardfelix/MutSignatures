import os
import trackhub

hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name="BRASS",
    genome="hg19",
    email="felix.racine-brassard@usherbrooke.ca")

track1 = trackhub.Track(
    name="HAP1_TP53_HU_150uM_C1",
    source='HAP1_TP53_HU_150uM_C1.BRASS.picardFILTERED.bam',
    tracktype='bam',
)

trackdb.add_tracks(track1)