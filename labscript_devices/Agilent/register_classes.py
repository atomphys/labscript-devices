
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'Agilent',
    BLACS_tab='labscript_devices.Agilent.blacs_tab.AgilentTab',
    runviewer_parser=None
)
