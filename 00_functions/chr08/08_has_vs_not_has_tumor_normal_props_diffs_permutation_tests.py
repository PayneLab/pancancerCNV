import cnvutils

cnvutils.props_permutation_test(
    n=10000,
    sources=["cptac", "gistic"],
    levels=["gene"],
    chromosomes_events={
        8: {
            "p": ["loss"],
            "q": ["gain"],
        },
    },
)