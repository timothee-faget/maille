//! Encodage PNG minimal : RGBA 8 bits, sans compression (blocs « stockés » de zlib). Les fichiers
//! sont gros mais valides, et ça évite une dépendance pour des captures d'écran.

/// PNG d'une image `width` × `height` dont `rgba` donne les pixels, ligne par ligne du haut vers le bas.
pub fn encode(width: u32, height: u32, rgba: &[u8]) -> Vec<u8> {
    assert_eq!(rgba.len(), (width * height * 4) as usize);
    let mut raw = Vec::with_capacity(rgba.len() + height as usize);
    for row in rgba.chunks((width * 4) as usize) {
        raw.push(0); // filtre « aucun »
        raw.extend_from_slice(row);
    }
    let mut zlib = vec![0x78, 0x01];
    let mut blocks = raw.chunks(65535).peekable();
    if blocks.peek().is_none() {
        zlib.extend_from_slice(&[1, 0, 0, 0xff, 0xff]);
    }
    while let Some(block) = blocks.next() {
        zlib.push(u8::from(blocks.peek().is_none()));
        zlib.extend_from_slice(&(block.len() as u16).to_le_bytes());
        zlib.extend_from_slice(&(!(block.len() as u16)).to_le_bytes());
        zlib.extend_from_slice(block);
    }
    zlib.extend_from_slice(&adler32(&raw).to_be_bytes());

    let mut png = b"\x89PNG\r\n\x1a\n".to_vec();
    let mut header = Vec::new();
    header.extend_from_slice(&width.to_be_bytes());
    header.extend_from_slice(&height.to_be_bytes());
    header.extend_from_slice(&[8, 6, 0, 0, 0]); // 8 bits, RGBA, sans entrelacement
    chunk(&mut png, b"IHDR", &header);
    chunk(&mut png, b"IDAT", &zlib);
    chunk(&mut png, b"IEND", &[]);
    png
}

fn chunk(png: &mut Vec<u8>, kind: &[u8; 4], data: &[u8]) {
    png.extend_from_slice(&(data.len() as u32).to_be_bytes());
    let start = png.len();
    png.extend_from_slice(kind);
    png.extend_from_slice(data);
    let crc = crc32(&png[start..]);
    png.extend_from_slice(&crc.to_be_bytes());
}

fn crc32(data: &[u8]) -> u32 {
    let mut crc = !0_u32;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xEDB8_8320 & (!(crc & 1)).wrapping_add(1));
        }
    }
    !crc
}

fn adler32(data: &[u8]) -> u32 {
    let (mut a, mut b) = (1_u32, 0_u32);
    for &byte in data {
        a = (a + u32::from(byte)) % 65521;
        b = (b + a) % 65521;
    }
    (b << 16) | a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checksums_match_their_reference_values() {
        assert_eq!(crc32(b"IEND"), 0xAE42_6082);
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(adler32(b"Wikipedia"), 0x11E6_0398);
    }

    /// Relit le PNG écrit : signature, IHDR, blocs stockés de zlib, somme de contrôle, IEND.
    #[test]
    fn the_image_can_be_read_back_pixel_for_pixel() {
        // plus de 65 535 octets de pixels : plusieurs blocs stockés
        let (w, h) = (200_u32, 100_u32);
        let rgba: Vec<u8> = (0..w * h * 4).map(|i| (i * 7 % 251) as u8).collect();
        let png = encode(w, h, &rgba);
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        let mut pos = 8;
        let mut idat = Vec::new();
        let mut kinds = Vec::new();
        while pos < png.len() {
            let len = u32::from_be_bytes(png[pos..pos + 4].try_into().unwrap()) as usize;
            let kind = &png[pos + 4..pos + 8];
            let data = &png[pos + 8..pos + 8 + len];
            let crc = u32::from_be_bytes(png[pos + 8 + len..pos + 12 + len].try_into().unwrap());
            assert_eq!(crc, crc32(&png[pos + 4..pos + 8 + len]), "CRC de {kind:?}");
            if kind == b"IHDR" {
                assert_eq!(u32::from_be_bytes(data[..4].try_into().unwrap()), w);
                assert_eq!(u32::from_be_bytes(data[4..8].try_into().unwrap()), h);
            }
            if kind == b"IDAT" {
                idat.extend_from_slice(data);
            }
            kinds.push(String::from_utf8(kind.to_vec()).unwrap());
            pos += 12 + len;
        }
        assert_eq!(kinds, ["IHDR", "IDAT", "IEND"]);
        // zlib : en-tête de 2 octets, blocs stockés, Adler-32
        assert_eq!(&idat[..2], &[0x78, 0x01]);
        let (mut raw, mut at, mut last) = (Vec::new(), 2, false);
        while !last {
            last = idat[at] == 1;
            let len = u16::from_le_bytes([idat[at + 1], idat[at + 2]]) as usize;
            assert_eq!(
                u16::from_le_bytes([idat[at + 3], idat[at + 4]]),
                !(len as u16)
            );
            raw.extend_from_slice(&idat[at + 5..at + 5 + len]);
            at += 5 + len;
        }
        assert_eq!(
            u32::from_be_bytes(idat[at..at + 4].try_into().unwrap()),
            adler32(&raw)
        );
        let rows: Vec<&[u8]> = raw.chunks((w * 4 + 1) as usize).collect();
        assert_eq!(rows.len(), h as usize);
        for (row, expected) in rows.iter().zip(rgba.chunks((w * 4) as usize)) {
            assert_eq!(row[0], 0);
            assert_eq!(&row[1..], expected);
        }
    }
}
