/* datauri

Create base64encoded `data` URIs from binary data.

Copyright (c) 2009 Sam Angove <sam [a inna circle] rephrase [period] net>

License: MIT

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/*
Base64 encode binary data for ASCII transport

See: http://en.wikipedia.org/wiki/Base64
*/
var base64encode = function(){
    // This is a non-standard extension available in Mozilla
    // and possibly other browsers.
    if (typeof window.btoa != "undefined")
        return window.btoa;

    /* JS fallback based on public domain code from Tyler Akins:
        http://rumkin.com/tools/compression/base64.php */
    var chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";
    function _btoa(data) {
        var chr1, chr2, chr3, enc1, enc2, enc3,
            i=0, length=data.length, output="";
        while (i < length) {
            // Convert 3 bytes of data into 4 6-bit chunks
            chr1 = data.charCodeAt(i++);
            chr2 = data.charCodeAt(i++);
            chr3 = data.charCodeAt(i++);

            enc1 = chr1 >> 2;                       // reduce byte to 6 bits
            enc2 = ((chr1 & 3) << 4) | (chr2 >> 4); // last 2 bits of chr1 + first 4 of chr2
            enc3 = ((chr2 & 15) << 2) | (chr3 >> 6);// last 4 bits of chr2 + first 2 of chr3
            enc4 = chr3 & 63;                       // last 6 bits

            if (isNaN(chr2)) enc3 = enc4 = 64;      // pad with zeroes if necessary
            else if (isNaN(chr3)) enc4 = 64;

            output += chars.charAt(enc1) + chars.charAt(enc2) + chars.charAt(enc3) + chars.charAt(enc4);
        }
        return output;
    }
    return _btoa;
}();

/*
Convert binary data to a `data` URI.

    See: http://en.wikipedia.org/wiki/Data_URI_scheme
*/
function datauri(content_type, data) {
    return "data:" + content_type
                   + ";base64,"
                   + encodeURIComponent(base64encode(data));
}