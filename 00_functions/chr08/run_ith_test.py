import cnvutils
import numpy as np
import os

ENTROPY = 294341852078515767793057943394766554310

# Get the index of this iteration
i = os.getenv("SLURM_ARRAY_TASK_ID")

# Generate the seed sequence for this iteration
sseq = np.random.SeedSequence(entropy=ENTROPY, spawn_key=(i,))

print(sseq.state)

# Create an RNG from that seed sequence
rng = np.random.default_rng(sseq)

# Get the results for this permutation
results = cnvutils.permute_props(
    sources=["cptac", "gistic"],
    levels=["gene"],
    chromosomes_events={
        8: {
            "p": ["loss"],
            "q": ["gain"],
        },
    },
    rng=rng,
)

# Save the results to a temporary output directory
results.to_csv(
    os.path.join("tmp", f"results_{i:0>6}"),
    sep="\t",
)
