# Spam Align

## Notes
The algorithm for checking and updating consistency bound is based on [GABIOS-LIB](gobics.de/burkhard/papers/jobim.pdf).  


## Pattern generation

A better way of generating patterns is needed. The current function with the signature
```rust
pub fn generate_random_patterns(
    count: usize,
    min_len: usize,
    max_len: usize,
    min_weight: usize,
    max_weight: usize,
) -> Vec<Pattern> {}
```

is not flexible enough. (One current limitation is that the min_len cannot be less than the max_weight).
A better approach might be to introduce a trait that which gets implemented by different pattern factories.  
Those factories could be more or less complicated and can then just be plugged into the alignment function in order to generate the patterns.  
A better alternative to the above function might be a factory that takes a count, a min lenght, a max length and two distributions, one for the distribution of the pattern lengths (should there be more short patterns than long ones?) and one for the occurrence of match positions in those patterns (should they be centered, evenly distributed, etc.)


Other Idea: Maybe some appropriate distribution could also be used to calculate the number of match positions depending on the pattern size. For smaller patterns it might be better to have less dont care positions than bigger patterns.

## Todos


- [ ] Rework pattern generation
    - maybe Rasbhari can be used, but i feel like a simpler solution would suffice for the moment
        - maybe save this for the thesis
    - see above Pattern generation section
- [ ] Rework generation of diagonals
    - at the moment this takes a considerable amount of time, this should be reduced by AT LEAST one order of magnitude
    - also diagonal generation should be parallelized (benchmarks needed!)
- [ ] provide writer for alignment results
    - which format? fasta would be simplest but sucks
- [ ] optimize updating trans bounds by utilizing eq classes like described in the linked paper at the top
- [ ] evaluate alignments against bb by computing appropriate scores
    - research needed about what is needed
- [ ] start to write this shit up! it will take longer than expected