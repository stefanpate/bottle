# from equilibrator_api import ComponentContribution, Q_
# cc = ComponentContribution()

#################################################################

from equilibrator_assets.local_compound_cache import LocalCompoundCache
lc = LocalCompoundCache()
lc.generate_local_cache_from_default_zenodo("compounds.sqlite")