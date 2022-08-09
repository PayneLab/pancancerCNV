import cnvutils
import numpy as np
import os
import sys

ENTROPY = 294341852078515767793057943394766554310

# Get the index of this iteration
i = int(sys.argv[1])

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
save_dir = f"tmp_raw{i // 1000:0>3}"
os.makedirs(save_dir, exist_ok=True)

results.to_csv(
    os.path.join(save_dir, f"results_{i:0>6}.tsv"),
    sep="\t",
)
