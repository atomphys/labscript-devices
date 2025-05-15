
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'data_logger',
    BLACS_tab='labscript_devices.data_logger.blacs_tab.data_loggerTab',
    runviewer_parser=None
)
