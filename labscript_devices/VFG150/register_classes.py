
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'VFG150',
    BLACS_tab='labscript_devices.VFG150.blacs_tab.VFG150Tab',
    runviewer_parser=None
)
