
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'Keysight_Scope',
    BLACS_tab='labscript_devices.Keysight_Scope.blacs_tab.Keysight_ScopeTab',
    runviewer_parser=None
)
