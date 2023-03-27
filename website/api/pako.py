# https://stackoverflow.com/questions/46351275/using-pako-deflate-with-python

from urllib.parse import quote, unquote
import base64
import zlib

class Js:
    @staticmethod
    def encode_uri_component(data):
        return quote(data, safe='~()*!.\'')

    @staticmethod
    def decode_uri_component(data):
        return unquote(data)

    @staticmethod
    def string_to_byte(data):
        return bytes(data, 'iso-8859-1')

    @staticmethod
    def bytes_to_string(data):
        return data.decode('iso-8859-1')

    @staticmethod
    def btoa(data):
        return base64.b64encode(data)

    @staticmethod
    def atob(data):
        return base64.b64decode(data)

class Pako:
    @staticmethod
    def deflate(data):
        compress  = zlib.compressobj(zlib.Z_DEFAULT_COMPRESSION, zlib.DEFLATED, 15, 
            memLevel=8, strategy=zlib.Z_DEFAULT_STRATEGY)
        compressed_data = compress.compress(Js.string_to_byte(Js.encode_uri_component(data)))
        compressed_data += compress.flush()
        return compressed_data

    @staticmethod
    def deflate_raw(data):
        compress = zlib.compressobj(
            zlib.Z_DEFAULT_COMPRESSION, zlib.DEFLATED, -15, memLevel=8,
            strategy=zlib.Z_DEFAULT_STRATEGY)
        compressed_data = compress.compress(Js.string_to_byte(Js.encode_uri_component(data)))
        compressed_data += compress.flush()
        return compressed_data
    
    @staticmethod
    def inflate(data):
        decompress = zlib.decompressobj(15)
        decompressed_data = decompress.decompress(data)
        decompressed_data += decompress.flush()
        return decompressed_data

    @staticmethod
    def inflate_raw(data):
        decompress = zlib.decompressobj(-15)
        decompressed_data = decompress.decompress(data)
        decompressed_data += decompress.flush()
        return decompressed_data