#!/bin/env python

##########################################################################
#	Object Class for the Metrics Dictionary
##########################################################################
from statistics import mean, median
import numpy as np # type: ignore

## initialize all the metrics
## with all the default values 
## then I don't have to worry about missing values 

def safe_mean(values):
        """
        Return mean of array accounting for None values or empty arrays
        """
        if values:
            filtered_values = [v for v in values if v is not None]
            return mean(filtered_values) if filtered_values is not None else None
        return 0 

def safe_median(values):
    """
    Return median of array accounting for None values or empty arrays
    """
    if values:
        filtered_values = [v for v in values if v is not None]
        return median(filtered_values) if filtered_values is not None else None
    return 0

class MetricsDictionary:
    """
    Class to hold the Base and Window Metrics dictionary
    """
    def __init__(self, sample, cohort, index):
        self.metrics = {
            "vcf_index" : index,
            "Sample" : sample,
            "Cohort" : cohort,
            "tumor_other_bases_count" : 0,
            "tumor_depth" : None,
            "tumor_VAF" : None,
            "num_total_reads" : 0,
            "tumor_reads_filtered": 0,
            "FDeamC" : 0,
            "SOB" : 0,
            "trinucleotide_context" : "",
            "pentanucleotide_context": "",
            "window_gc_cont": 0,
            "window_seq_entropy": 0,
            "window_median_cov": 0 ,
            "window_cov_variance" : 0,
            "window_min_cov_ratio": 0,
            "window_max_cov_ratio": 0,
            "window_median_frag_len": 0,
            "window_dup_frac": 0,
            "window_multi_frac": 0,
            "window_improper_frac" : 0,
            "window_median_mapq" : 0,
            "window_read_filter_frac" : 0,
            "tumor_ref_base_quality_frac" : 0
        }

        for prefix in ["tumor_ref", "tumor_var"]:
            self.metrics.update({
                f'{prefix}_count' : 0,
                f'{prefix}_num_plus_strand': 0,
                f'{prefix}_num_minus_strand': 0,
                f'{prefix}_base_qualities': [],
                f'{prefix}_mapping_quality': [],
                f'{prefix}_avg_num_mismatches': [],
                f'{prefix}_avg_sum_mismatch_base_quals': [],
                f'{prefix}_read_frag_length': [],
                f'{prefix}_clipped_length': [],
                f'{prefix}_avg_pos_as_fraction': [],
                f'{prefix}_distances_to_3p_end': [],
                f'{prefix}_distances_to_5p_end': []
            })

    def increment_metric(self, metric_name):
        """
        Increment a metric by 1 (useful for tumor_other_bases_count, num_total_reads,
        {prefix}_num_minus_strand, {prefix}_num_plus_strand)
        """
        if metric_name not in self.metrics:
            self.metrics[metric_name] = 0
        self.metrics[metric_name] += 1

    def add_metric(self, metric_name, value):
        """
        Add a metric to the dictionary
        """
        if metric_name not in self.metrics:
            self.metrics[metric_name] = []
        self.metrics[metric_name].append(value)

    def set_metric(self, metric_name, value):
        self.metrics[metric_name] = value

    def get_metric(self, metric_name):
        """
        Get a metric from the dictionary
        """
        return self.metrics.get(metric_name, [])

    def get_all_metrics(self):
        """
        Get all metrics from the dictionary
        """
        return self.metrics
    
    def update_metrics(self, prefix):
        """
        Factored code for updating pooled dictionary of features for all variants (metrics). 
        Will append pooled dictionary with newly calculated features. 

        Args:
            prefix (str): var or ref prefix
        """    
        ## Need to make sure I don't pop things before I use them
        self.metrics.update({
            f'{prefix}_avg_base_quality': safe_mean(self.metrics.get(f'{prefix}_base_qualities', [])),
            f'{prefix}_med_base_quality': safe_median(self.metrics.pop(f'{prefix}_base_qualities', [])),
            f'{prefix}_avg_mapping_quality': safe_mean(self.metrics.get(f'{prefix}_mapping_quality', [])),
            f'{prefix}_med_mapping_quality': safe_median(self.metrics.pop(f'{prefix}_mapping_quality', [])),
            f'{prefix}_avg_num_mismatches_as_fraction': safe_mean(self.metrics.pop(f'{prefix}_avg_num_mismatches', [])),
            f'{prefix}_avg_sum_mismatch_qualities': safe_mean(self.metrics.pop(f'{prefix}_avg_sum_mismatch_base_quals', [])),
            f'{prefix}_med_frag_len': safe_median(self.metrics.pop(f'{prefix}_read_frag_length', [])),
            f'{prefix}_avg_clipped_length': safe_mean(self.metrics.pop(f'{prefix}_clipped_length', [])),
            f'{prefix}_avg_pos_as_fraction': safe_mean(self.metrics.pop(f'{prefix}_avg_pos_as_fraction', [])),
            f'{prefix}_avg_distance_to_effective_3p_end': safe_mean(self.metrics.get('{prefix}_distances_to_3p_end', [])),
            f'{prefix}_avg_distance_to_effective_5p_end': safe_mean(self.metrics.get(f'{prefix}_distances_to_5p_end', [])),
            f'{prefix}_med_distance_to_effective_3p_end': safe_median(self.metrics.pop(f'{prefix}_distances_to_3p_end', [])),
            f'{prefix}_med_distance_to_effective_5p_end': safe_median(self.metrics.pop(f'{prefix}_distances_to_5p_end', []))})
        
        self.metrics.update({
            f'{prefix}_normmed_distance_to_effective_3p_end': self.metrics.get(f'{prefix}_med_distance_to_effective_3p_end')\
                / self.metrics[f'{prefix}_med_frag_len'] if self.metrics[f'{prefix}_med_frag_len'] != 0 else 0,
            f'{prefix}_normmed_distance_to_effective_5p_end': self.metrics[f'{prefix}_med_distance_to_effective_5p_end']\
                / self.metrics[f'{prefix}_med_frag_len'] if self.metrics[f'{prefix}_med_frag_len'] != 0 else 0
        })

        
    def aggregate_base_metrics(self, REF, ALT):
        self.update_metrics('tumor_ref')
        self.update_metrics('tumor_var')

        #what to do about this because i reversed it originally...
        # I guess I could just negate the FDEAMC calculations

        tumor_var_num_plus_strand = self.metrics.get('tumor_var_num_plus_strand', 0)
        tumor_var_num_minus_strand = self.metrics.get('tumor_var_num_minus_strand', 0)
        
        if tumor_var_num_plus_strand + tumor_var_num_minus_strand > 0:
            if REF == "C" and ALT == "T":
                self.metrics['FDeamC'] = float((-1) * tumor_var_num_plus_strand /\
                                                (tumor_var_num_plus_strand + \
                                                tumor_var_num_minus_strand))
            elif REF == "G" and ALT == "A":
                self.metrics['FDeamC'] =  float((-1) * tumor_var_num_minus_strand /\
                        (tumor_var_num_plus_strand + \
                            tumor_var_num_minus_strand))
            
            self.metrics['SOB'] = float((tumor_var_num_plus_strand - tumor_var_num_minus_strand) \
                                / (tumor_var_num_plus_strand + tumor_var_num_minus_strand))
        
    
        ref_base_qualities = self.metrics.get('tumor_ref_med_base_quality', 0)
        var_base_qualities = self.metrics.get('tumor_var_med_base_quality', 0)

        self.metrics['tumor_ref_base_quality_frac'] = (var_base_qualities / ref_base_qualities) if (ref_base_qualities > 0) else 0
        
        if self.metrics['tumor_VAF'] is None:
            self.metrics['tumor_VAF'] = self.metrics['tumor_var_count'] / self.metrics['tumor_depth'] if self.metrics['tumor_depth'] > 0 else 0            

    def update_coverage_ratios(self, left, right):
        window_coverage_ratios = (float(left / self.metrics['window_median_cov']), float(right / self.metrics['window_median_cov']))
        self.metrics['window_min_cov_ratio'] = min(window_coverage_ratios)
        self.metrics['window_max_cov_ratio'] = max(window_coverage_ratios)
