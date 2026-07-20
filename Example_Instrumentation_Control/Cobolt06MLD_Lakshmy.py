#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 09:21:39 2025

@author: lakshmypriyaajayakumar
"""

import serial
from serial.tools import list_ports
from serial.serialutil import SerialException
import time
import sys
import re
import logging

logger = logging.getLogger(__name__)

class CoboltLaser:
    """Creates a laser object using either COM-port or serial number to connect to laser. \n Will automatically return proper subclass, if applicable"""

    def __init__(
        self,
        port = 'COM3'',
        serialnumber = None,
        baudrate: int = 115200,
    ):
        self.msg_timer = time.perf_counter()
        self.serialnumber = str(serialnumber)
        self.port = port
        self.modelnumber = None
        self.baudrate = baudrate
        self.address = None
        self.connect()

    def __repr__(self):
        try:
            return f'Serial number: {self.serialnumber}, Model number: {self.modelnumber}, Wavelength: {"{:.0f}".format(float(self.modelnumber[0:4]))} nm, Type: {self.__class__.__name__} Port: {self.port}'
        except:
            return f"Serial number: {self.serialnumber}, Model number: {self.modelnumber}, Port: {self.port}"

    def connect(self):
        """Connects the laser on using a specified COM-port (preferred) or serial number. Will throw exception if it cannot connect to specified port or find laser with given serial number.

        Raises:
            SerialException: serial port error
            RuntimeError: no laser found
        """

        if self.port != None:
            try:
                self.address = serial.Serial(self.port, self.baudrate, timeout=1)
            except Exception as err:
                self.address = None
                raise SerialException(f"{self.port} not accesible.") from err

        elif self.serialnumber != None:
            ports = [x for x in list_ports.comports() if "USB" in x.hwid]
            try:
                port=next([item for item in ports if self.serialnumber in item.serial_number])
                self.port=port.device
                self.address = serial.Serial(self.port, self.baudrate, timeout=1)
            except:
                for port in ports:
                    try:
                        self.address = serial.Serial(
                            port.device, baudrate=self.baudrate, timeout=1
                        )
                        sn = self.send_cmd("sn?")
                        self.address.close()
                        if sn == self.serialnumber:
                            self.port = port.device
                            self.address = serial.Serial(self.port, baudrate=self.baudrate)
                            break
                    except:
                        pass
            if self.port == None:
                raise RuntimeError("No laser found")
        if self.address != None:
            self._identify_()
        if self.__class__ == CoboltLaser:
            self._classify_()

    def _identify_(self):
        """Fetch Serial number and model number of laser. Will raise exception and close connection if not connected to a cobolt laser.

        Raises:
            RuntimeError: error identifying the laser model
        """
        try:
            firmware = self.send_cmd("gfv?")
            if "error" in firmware.lower():
                self.disconnect()
                raise RuntimeError("Not a Cobolt laser")
            self.serialnumber = self.send_cmd("sn?")
            if not "." in firmware:
                if "0" in self.serialnumber:
                    self.modelnumber = (
                        f"0{self.serialnumber.partition(str(0))[0]}-04-XX-XXXX-XXX"
                    )
                    self.serialnumber = self.serialnumber.partition("0")[2]
                    while self.serialnumber[0] == "0":
                        self.serialnumber = self.serialnumber[1:]
            else:
                self.modelnumber = self.send_cmd("glm?")
        except:
            self.disconnect()
            raise RuntimeError("Not a Cobolt laser")

    def _classify_(self):
        """Classifies the laser into probler subclass depending on laser type"""
        try:
            if re.search("-06-.*-(1\d{3})(|-C)$", self.modelnumber):
                self.__class__ = Cobolt06
            elif re.search("-06-0.*(|-C)$", self.modelnumber):
                self.__class__ = Cobolt06MLD
            elif re.search("-06-(5|9).*-(\d{3})(|-C)$", self.modelnumber):
                self.__class__ = Cobolt06DPL
        except:
            pass

    def is_connected(self):
        """Ask if laser is connected"""
        try:
            if self.address.is_open:
                try:
                    test = self.send_cmd("?")
                    if test == "OK":
                        return True
                    else:
                        return False
                except:
                    return False
            else:
                return False
        except:
            return False

    def disconnect(self):
        """Disconnect the laser"""
        if self.address != None:
            self.address.close()
            self.serialnumber = None
            self.modelnumber = None

    def turn_on(self):
        """Turn on the laser with the autostart sequence.The laser will await the TEC setpoints and pass a warm-up state"""
        logger.info("Turning on laser")
        return self.send_cmd(f"@cob1")

    def turn_off(self):
        """Turn off the laser"""
        logger.info("Turning off laser")
        return self.send_cmd(f"l0")

    def is_on(self):
        """Ask if laser is turned on"""
        answer = self.send_cmd(f"l?")
        if answer == "1":
            return True
        else:
            return False

    def interlock(self):
        """Returns: 0 if closed, 1 if open"""
        return self.send_cmd(f"ilk?")

    def get_fault(self):
        """Get laser fault"""
        return self.send_cmd("f?")

    def clear_fault(self):
        """Clear laser fault"""
        return self.send_cmd("cf")

    def get_mode(self):
        """Get operating mode"""
        mode = self.send_cmd("gam?")
        return mode

    def get_state(self):
        """Get autostart state"""
        state = self.send_cmd("gom?")
        return state

    def constant_current(self, current=None):
        """Enter constant current mode, current in mA"""
        if current != None:
            if not "-08-" in self.modelnumber or not "-06-" in self.modelnumber:
                self.send_cmd(f"slc {current/1000}")
            else:
                self.send_cmd(f"slc {current}")
            logger.info(f"Entering constant current mode with I = {current} mA")
        else:
            logger.info("Entering constant current mode")
        return self.send_cmd(f"ci")

    def set_current(self, current:float):
        """Set laser current in mA"""
        logger.info(f"Setting I = {current} mA")
        if not "-08-" in self.modelnumber or not "-06-" in self.modelnumber:
            current = current / 1000
        return self.send_cmd(f"slc {current}")

    def get_current(self):
        """Get laser current in mA"""
        return float(self.send_cmd(f"i?"))

    def get_current_setpoint(self):
        """Get laser current setpoint in mA"""
        return float(self.send_cmd(f"glc?"))

    def constant_power(self, power:float|None=None):
        """Enter constant power mode, power in mW"""
        if power != None:
            self.send_cmd(f"p {float(power)/1000}")
            logger.info(f"Entering constant power mode with P = {power} mW")
        else:
            logger.info("Entering constant power mode")
        return self.send_cmd(f"cp")

    def set_power(self, power:float):
        """Set laser power in mW"""
        logger.info(f"Setting P = {power} mW")
        return self.send_cmd(f"p {float(power)/1000}")

    def get_power(self):
        """Get laser power in mW"""
        return float(self.send_cmd(f"pa?")) * 1000

    def get_power_setpoint(self):
        """Get laser power setpoint in mW"""
        return float(self.send_cmd(f"p?")) * 1000

    def get_ophours(self):
        """Get laser operational hours"""
        return self.send_cmd(f"hrs?")


    def send_cmd(self, message, timeout= 1):
        """Sends a message to the laser and awaits response until timeout (in s).

        Returns:
            The response received from the laser as string

        Raises:
            RuntimeError: sending the message failed
        """


        if timeout:
            self.address.timeout = timeout
        message += "\r"
        while time.perf_counter()-self.msg_timer<0.100: #to prevent issues with sending commands too rapidly 
            continue
        try:
            utf8_msg = message.encode()
            self.address.write(utf8_msg)
            logger.debug(f"sent laser [{self}] message [{utf8_msg}]")
        except Exception as e:
            raise RuntimeError("Error: write failed") from e

        try:
            received_string = self.address.readline().decode().rstrip()
            self.msg_timer = time.perf_counter()
            if len(received_string) < 1:  # if empty response raise syntax error
                logger.error(f"No responce recieved for {message}")
                raise SerialException
        except serial.SerialException:

            raise RuntimeError(f"Syntax Error: No response on {message}")
        else:
            #print(message.replace("\r",""),received_string)
            logger.debug(f"received from laser [{self}] message [{received_string}]")
        return received_string

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.turn_off()
        self.disconnect()
        
        
     """ For Cobolt 06-MLD 375 nm    """
        
        

class Cobolt06MLD(CoboltLaser):
    """For lasers of type 06-MLD"""

    def __init__(self, port=None, serialnumber=None):
        super().__init__(port, serialnumber)

    def get_mode(self):
        """Get operating mode"""
        modes = {
            "0": "0 - Constant Current",
            "1": "1 - Constant Power",
            "2": "2 - Modulation Mode",
        }
        mode = self.send_cmd("gam?")
        return modes.get(mode, mode)

    def get_state(self):
        """Get autostart state"""
        states = {
            "0": "0 - Off",
            "1": "1 - Waiting for key",
            "2": "2 - Continuous",
            "3": "3 - On/Off Modulation",
            "4": "4 - Modulation",
            "5": "5 - Fault",
            "6": "6 - Aborted",
        }
        state = self.send_cmd("gom?")
        return states.get(state, state)

    def modulation_mode(self, power:float=None):
        """Enter modulation mode.

        Args:
            power: modulation power (mW)
        """
        logger.info(f"Entering modulation mode")
        if power != None:
            self.send_cmd(f"slmp {power}")
        return self.send_cmd("em")

    def digital_modulation(self, enable:int):
        """Enable digital modulation mode by enable=1, turn off by enable=0"""
        return self.send_cmd(f"sdmes {enable}")

    def analog_modulation(self, enable:int):
        """Enable analog modulation mode by enable=1, turn off by enable=0"""
        return self.send_cmd(f"sames {enable}")

    def on_off_modulation(self, enable:int):
        """Enable On/Off modulation mode by enable=1, turn off by enable=0"""
        if enable == 1:
            return self.send_cmd("eoom")
        elif enable == 0:
            return self.send_cmd("xoom")

    def get_modulation_state(self):
        """Get the laser modulation settings as [analog, digital]"""
        dm = self.send_cmd("gdmes?")
        am = self.send_cmd("games?")
        return [am, dm]

    def set_modulation_power(self, power:float):
        """Set the modulation power in mW"""
        logger.info(f"Setting modulation power = {power} mW")
        return self.send_cmd(f"slmp {power}")

    def get_modulation_power(self):
        """Get the modulation power setpoint in mW"""
        return float(self.send_cmd("glmp?"))

    def set_analog_impedance(self, arg:int):
        """Set the impedance of the analog modulation.

        Args:
            arg: 0 for HighZ, 1 for 50 Ohm.
        """
        return self.send_cmd(f"salis {arg}")

    def get_analog_impedance(self):
        """Get the impedance of the analog modulation \n
        return: 0 for HighZ and 1 for 50 Ohm"""
        return self.send_cmd("galis?")
