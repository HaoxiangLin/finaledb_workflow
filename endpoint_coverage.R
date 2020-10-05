# The MIT License (MIT)
# 
# Copyright (c) 2011-2020 Cincinnati Children's Hospital Medical Center
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
# Purpose: for each genome coordinate, calculate the number of fragments with 
# endpoints at this locus.
#
# Usage: endpoints_coverage.R fragment.tsv.gz endpoints_coverage_start.bed.gz endpoints_coverage_end.bed.gz
#
# Input: a gzipped file of tab-separated values format (TSV), with the following
# columns: chromosome name, fragment start position, fragment end position, MAPQ,
# and strand specification.

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
frags_file <- args[1]
output_file_start <- args[2]
output_file_end <- args[3]

frags <-
  read_tsv(
    frags_file,
    col_names = c("chr", "start", "end", "mapq", "strand"),
    col_types = cols(
      col_character(),
      col_integer(),
      col_integer(),
      col_integer(),
      col_character()
    )
  ) %>%
  mutate(chr = factor(chr, levels = c(1:22, "X", "Y")))

frags_long <- frags %>% select(chr, start, end) %>% pivot_longer(c(start, end), names_to = "type")

endpoints_coverage_start <- frags_long %>% 
  filter(type == "start") %>% 
  rename(start = value) %>% 
  group_by(chr, start) %>% 
  summarize(coverage = length(start), .groups = "drop") %>%
  mutate(end = start + 1) %>%
  arrange(chr, start) %>%
  select(chr, start, end, coverage)


endpoints_coverage_end <- frags_long %>% 
  filter(type == "end") %>% 
  rename(end = value) %>% 
  group_by(chr, end) %>% 
  summarize(coverage = length(end), .groups = "drop") %>%
  mutate(start = end - 1) %>%
  arrange(chr, start) %>%
  select(chr, start, end, coverage)

write_tsv(endpoints_coverage_start, output_file_start, col_names = FALSE)
write_tsv(endpoints_coverage_end, output_file_end, col_names = FALSE)
