
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'iq_modulator',
    BLACS_tab='labscript_devices.iq_modulator.blacs_tab.iq_modulatorTab',
    runviewer_parser=None
)
