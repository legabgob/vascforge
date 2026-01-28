from vascx.fundus.features.bifurcation_counts import BifurcationCount
from vascx.fundus.features.caliber import Caliber
from vascx.fundus.features.cre import CRE
from vascx.fundus.features.temporal_angles import TemporalAngle
from vascx.fundus.features.tortuosity import LengthMeasure, Tortuosity
from vascx.fundus.features.vascular_densities import VascularDensity
from vascx.shared.aggregators import median, median_std
from vascx.shared.bl_aggregators import root, middle, leaves
from vascx.shared.features import FeatureSet
from vascx.fundus.features.regression_slopes import RegressionSlope

bergmann_features = FeatureSet(
    "bergmann",
    {
        # "ta": TemporalAngle(),
        # "cre": CRE(),
        # "vd": VascularDensity(),
        # "diam": Caliber(aggregator=median_std),
        #"tort": Tortuosity(length_measure=LengthMeasure.Skeleton, aggregator=median),
        # "bif": BifurcationCount(),
        #'bif_root': BifurcationCount(bl_aggregator=root),
        #'bif_middle': BifurcationCount(bl_aggregator=middle),
        #'bif_leaves': BifurcationCount(bl_aggregator=leaves),
        'tort_root': Tortuosity(length_measure=LengthMeasure.Splines, aggregator=median, bl_aggregator=root),
        'tort_middle': Tortuosity(length_measure=LengthMeasure.Splines, aggregator=median, bl_aggregator=middle),
        'tort_leaves': Tortuosity(length_measure=LengthMeasure.Splines, aggregator=median, bl_aggregator=leaves),
        #'diam_root': Caliber(aggregator=median_std, bl_aggregator=root),
        #'diam_middle': Caliber(aggregator=median_std, bl_aggregator=middle),
        #'diam_leaves': Caliber(aggregator=median_std, bl_aggregator=leaves),
        'regr_slope_diam': RegressionSlope(feature='diameter', levels_type='parts', r_squared=False),
        'regr_slope_tort': RegressionSlope(feature='tortuosity', levels_type='parts', r_squared=False),
        'regr_slope_bif': RegressionSlope(feature='bifurcation', levels_type='parts', r_squared=False),
    },
)
