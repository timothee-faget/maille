//! Lecture du format texte d'un STEP (ISO 10303-21, « Part 21 ») : en-tête, puis instances
//! numérotées `#12 = NOM(paramètres);` de la section DATA, simples ou complexes
//! (`#12 = (NOM_A(…) NOM_B(…));`). Les paramètres sont gardés tels quels (entiers, réels, chaînes,
//! énumérations, références, listes, valeurs typées) ; leur sens est l'affaire de la topologie.

use std::collections::BTreeMap;

use super::error;
use crate::IoError;

/// Contenu d'un fichier Part 21.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Part21 {
    /// Schémas de `FILE_SCHEMA`, par exemple `AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }` (AP214).
    pub schemas: Vec<String>,
    /// Instances de la section DATA, par numéro.
    pub instances: BTreeMap<u64, Instance>,
}

/// Instance `#n` : un enregistrement (instance simple) ou plusieurs (instance complexe, dans
/// l'ordre du fichier).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Instance {
    /// Ligne de `#n` dans le fichier, pour les messages d'erreur.
    pub line: usize,
    pub records: Vec<Record>,
}

/// `NOM(paramètres)`.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Record {
    pub name: String,
    pub params: Vec<Value>,
}

/// Paramètre d'un enregistrement.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Value {
    Integer(i64),
    Real(f64),
    String(String),
    /// `.NOM.`, sans les points (`T`, `F`, `UNSPECIFIED`…).
    Enum(String),
    /// `#n`.
    Ref(u64),
    /// `$` : valeur absente.
    Null,
    /// `*` : valeur dérivée.
    Derived,
    List(Vec<Value>),
    /// `NOM(valeur)`, par exemple `LENGTH_MEASURE(2.5)`.
    Typed(String, Vec<Value>),
    /// `"…"` : chaîne de bits en hexadécimal, gardée telle quelle.
    Binary(String),
}

impl Part21 {
    /// Instance `#id`, ou une erreur qui dit qu'elle manque (citée à la ligne `line`).
    pub fn get(&self, id: u64, line: usize) -> Result<&Instance, IoError> {
        self.instances.get(&id).ok_or_else(|| {
            error(format!(
                "ligne {line} : l'instance #{id} est citée mais n'existe pas dans la section DATA"
            ))
        })
    }
}

impl Instance {
    /// Enregistrement nommé `name` (l'unique d'une instance simple, l'un de ceux d'une complexe).
    pub fn record(&self, name: &str) -> Option<&Record> {
        self.records.iter().find(|r| r.name == name)
    }
}

/// Lit un fichier Part 21.
pub(super) fn parse(text: &str) -> Result<Part21, IoError> {
    let trimmed = text.trim_start();
    if trimmed.starts_with('<') {
        return Err(error(
            "STEP en XML (ISO 10303-28) : Maille lit le format texte ISO 10303-21, ré-exporter en .step ou .stp"
                .to_string(),
        ));
    }
    if trimmed.starts_with("PK") {
        return Err(error(
            "STEP compressé (.stpZ ou archive) : le décompresser d'abord".to_string(),
        ));
    }
    let mut p = Parser {
        lexer: Lexer::new(text),
        peeked: None,
    };
    p.keyword(
        "ISO-10303-21",
        "le fichier ne commence pas par « ISO-10303-21; » : ce n'est pas un STEP",
    )?;
    p.punct(b';')?;
    p.keyword("HEADER", "section HEADER attendue")?;
    p.punct(b';')?;
    let mut schemas = Vec::new();
    loop {
        let (line, name) = p.name()?;
        if name == "ENDSEC" {
            p.punct(b';')?;
            break;
        }
        let params = p.params(line)?;
        p.punct(b';')?;
        if name == "FILE_SCHEMA" {
            if let Some(Value::List(list)) = params.first() {
                schemas = list
                    .iter()
                    .filter_map(|v| match v {
                        Value::String(s) => Some(s.clone()),
                        _ => None,
                    })
                    .collect();
            }
        }
    }

    let mut instances = BTreeMap::new();
    loop {
        let (line, name) = p.name()?;
        match name.as_str() {
            "END-ISO-10303-21" => {
                p.punct(b';')?;
                break;
            }
            "DATA" => {
                // DATA ou DATA(nom, schéma) d'un fichier à plusieurs sections
                if p.peek_punct(b'(')? {
                    p.params(line)?;
                }
                p.punct(b';')?;
                p.data(&mut instances)?;
            }
            // sections facultatives qu'on ne lit pas
            "ANCHOR" | "REFERENCE" | "SIGNATURE" => p.skip_section()?,
            _ => {
                return Err(error(format!(
                    "ligne {line} : section DATA ou fin « END-ISO-10303-21; » attendue, trouvé « {name} »"
                )))
            }
        }
    }
    if instances.is_empty() {
        return Err(error(
            "aucune instance dans la section DATA : le fichier ne décrit aucune pièce".to_string(),
        ));
    }
    Ok(Part21 { schemas, instances })
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// Mot-clé : `CARTESIAN_POINT`, `ISO-10303-21`, `ENDSEC`…
    Keyword(String),
    Instance(u64),
    Integer(i64),
    Real(f64),
    String(String),
    Enum(String),
    Binary(String),
    Punct(u8),
    End,
}

struct Lexer<'a> {
    bytes: &'a [u8],
    pos: usize,
    line: usize,
}

impl<'a> Lexer<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            bytes: text.as_bytes(),
            pos: 0,
            line: 1,
        }
    }

    fn at(&self, k: usize) -> Option<u8> {
        self.bytes.get(self.pos + k).copied()
    }

    fn bump(&mut self) -> Option<u8> {
        let c = self.at(0)?;
        self.pos += 1;
        if c == b'\n' {
            self.line += 1;
        }
        Some(c)
    }

    /// Saute les blancs et les commentaires `/* … */`.
    fn skip(&mut self) -> Result<(), IoError> {
        loop {
            match self.at(0) {
                Some(c) if c.is_ascii_whitespace() => {
                    self.bump();
                }
                Some(b'/') if self.at(1) == Some(b'*') => {
                    let line = self.line;
                    self.pos += 2;
                    loop {
                        match self.bump() {
                            None => {
                                return Err(error(format!(
                                    "ligne {line} : commentaire « /* » jamais refermé"
                                )))
                            }
                            Some(b'*') if self.at(0) == Some(b'/') => {
                                self.pos += 1;
                                break;
                            }
                            Some(_) => {}
                        }
                    }
                }
                _ => return Ok(()),
            }
        }
    }

    /// Prochain lexème et sa ligne.
    fn next(&mut self) -> Result<(usize, Token), IoError> {
        self.skip()?;
        let line = self.line;
        let Some(c) = self.at(0) else {
            return Ok((line, Token::End));
        };
        let token = match c {
            b'(' | b')' | b',' | b'=' | b';' | b'$' | b'*' => {
                self.pos += 1;
                Token::Punct(c)
            }
            b'#' => {
                self.pos += 1;
                let digits = self.take(|c| c.is_ascii_digit());
                Token::Instance(digits.parse().map_err(|_| {
                    error(format!(
                        "ligne {line} : numéro d'instance attendu après « # »"
                    ))
                })?)
            }
            b'\'' => Token::String(self.string(line)?),
            b'"' => {
                self.pos += 1;
                let bits = self.take(|c| c.is_ascii_hexdigit());
                if self.bump() != Some(b'"') {
                    return Err(error(format!(
                        "ligne {line} : chaîne binaire « \" » mal fermée"
                    )));
                }
                Token::Binary(bits)
            }
            b'.' => {
                self.pos += 1;
                let name = self.take(|c| c.is_ascii_alphanumeric() || c == b'_');
                if name.is_empty() || self.bump() != Some(b'.') {
                    return Err(error(format!(
                        "ligne {line} : énumération « .NOM. » mal formée"
                    )));
                }
                Token::Enum(name)
            }
            b'+' | b'-' | b'0'..=b'9' => self.number(line)?,
            c if c.is_ascii_alphabetic() || c == b'!' => Token::Keyword(
                self.take(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-' || c == b'!'),
            ),
            c => {
                return Err(error(format!(
                    "ligne {line} : caractère inattendu « {} »",
                    c as char
                )))
            }
        };
        Ok((line, token))
    }

    fn take(&mut self, keep: impl Fn(u8) -> bool) -> String {
        let start = self.pos;
        while self.at(0).is_some_and(&keep) {
            self.pos += 1;
        }
        String::from_utf8_lossy(&self.bytes[start..self.pos]).into_owned()
    }

    /// Entier ou réel (un réel a un point ou un exposant : `1.`, `-2.5E-03`).
    fn number(&mut self, line: usize) -> Result<Token, IoError> {
        let start = self.pos;
        if matches!(self.at(0), Some(b'+' | b'-')) {
            self.pos += 1;
        }
        let mut real = false;
        while let Some(c) = self.at(0) {
            match c {
                b'0'..=b'9' => {}
                b'.' => real = true,
                b'E' | b'e' => {
                    real = true;
                    if matches!(self.at(1), Some(b'+' | b'-')) {
                        self.pos += 1;
                    }
                }
                _ => break,
            }
            self.pos += 1;
        }
        let text = std::str::from_utf8(&self.bytes[start..self.pos]).unwrap_or("");
        let bad = || error(format!("ligne {line} : nombre mal formé « {text} »"));
        if real {
            text.parse().map(Token::Real).map_err(|_| bad())
        } else {
            text.parse().map(Token::Integer).map_err(|_| bad())
        }
    }

    /// Chaîne `'…'` : `''` est une apostrophe ; directives `\\`, `\S\c` (ISO 8859-1 haut),
    /// `\X\hh`, `\X2\hhhh…\X0\` (UTF-16) et `\X4\hhhhhhhh…\X0\` ; `\Pc\` (page de code) ignorée. Les
    /// fins de ligne dans une chaîne ne comptent pas (une chaîne longue est coupée à l'écriture).
    fn string(&mut self, line: usize) -> Result<String, IoError> {
        let unclosed = || error(format!("ligne {line} : chaîne « ' » jamais refermée"));
        let bad = |what: &str| {
            error(format!(
                "ligne {line} : directive {what} mal formée dans une chaîne"
            ))
        };
        self.pos += 1;
        let mut raw = Vec::new();
        loop {
            match self.bump().ok_or_else(unclosed)? {
                b'\'' if self.at(0) == Some(b'\'') => {
                    self.pos += 1;
                    raw.push(b'\'');
                }
                b'\'' => break,
                b'\n' | b'\r' => {}
                c => raw.push(c),
            }
        }
        let text = String::from_utf8_lossy(&raw).into_owned();
        let mut out = String::new();
        let mut rest = text.as_str();
        while let Some(k) = rest.find('\\') {
            out.push_str(&rest[..k]);
            rest = &rest[k..];
            let hex = |s: &str| u32::from_str_radix(s, 16).ok();
            if let Some(r) = rest.strip_prefix("\\\\") {
                out.push('\\');
                rest = r;
            } else if let Some(r) = rest.strip_prefix("\\S\\") {
                let c = r.chars().next().ok_or_else(|| bad("\\S\\"))?;
                out.push(char::from_u32(c as u32 + 128).ok_or_else(|| bad("\\S\\"))?);
                rest = &r[c.len_utf8()..];
            } else if let Some(r) = rest.strip_prefix("\\X2\\").or(rest.strip_prefix("\\X4\\")) {
                let width = if rest.starts_with("\\X2\\") { 4 } else { 8 };
                let end = r.find("\\X0\\").ok_or_else(|| bad("\\X2\\ ou \\X4\\"))?;
                let digits = &r[..end];
                if digits.len() % width != 0 {
                    return Err(bad("\\X2\\ ou \\X4\\"));
                }
                let units: Vec<u32> = (0..digits.len() / width)
                    .map(|i| hex(&digits[i * width..(i + 1) * width]))
                    .collect::<Option<_>>()
                    .ok_or_else(|| bad("\\X2\\ ou \\X4\\"))?;
                if width == 4 {
                    let units: Vec<u16> = units.iter().map(|&u| u as u16).collect();
                    out.push_str(&String::from_utf16(&units).map_err(|_| bad("\\X2\\"))?);
                } else {
                    for u in units {
                        out.push(char::from_u32(u).ok_or_else(|| bad("\\X4\\"))?);
                    }
                }
                rest = &r[end + 4..];
            } else if let Some(r) = rest.strip_prefix("\\X\\") {
                let code = r.get(..2).and_then(hex).ok_or_else(|| bad("\\X\\"))?;
                out.push(char::from_u32(code).ok_or_else(|| bad("\\X\\"))?);
                rest = &r[2..];
            } else if rest.starts_with("\\P") && rest.get(3..4) == Some("\\") {
                rest = &rest[4..];
            } else {
                out.push('\\');
                rest = &rest[1..];
            }
        }
        out.push_str(rest);
        Ok(out)
    }
}

struct Parser<'a> {
    lexer: Lexer<'a>,
    peeked: Option<(usize, Token)>,
}

impl Parser<'_> {
    fn next(&mut self) -> Result<(usize, Token), IoError> {
        match self.peeked.take() {
            Some(t) => Ok(t),
            None => self.lexer.next(),
        }
    }

    fn peek(&mut self) -> Result<&(usize, Token), IoError> {
        if self.peeked.is_none() {
            self.peeked = Some(self.lexer.next()?);
        }
        Ok(self.peeked.as_ref().unwrap_or(&(0, Token::End)))
    }

    fn peek_punct(&mut self, c: u8) -> Result<bool, IoError> {
        Ok(self.peek()?.1 == Token::Punct(c))
    }

    fn unexpected(line: usize, expected: &str, found: &Token) -> IoError {
        let found = match found {
            Token::Keyword(k) => format!("« {k} »"),
            Token::Instance(n) => format!("« #{n} »"),
            Token::Punct(c) => format!("« {} »", *c as char),
            Token::End => "la fin du fichier".to_string(),
            _ => "une valeur".to_string(),
        };
        error(format!("ligne {line} : {expected} attendu, trouvé {found}"))
    }

    fn punct(&mut self, c: u8) -> Result<(), IoError> {
        match self.next()? {
            (_, Token::Punct(d)) if d == c => Ok(()),
            (line, t) => Err(Self::unexpected(line, &format!("« {} »", c as char), &t)),
        }
    }

    fn name(&mut self) -> Result<(usize, String), IoError> {
        match self.next()? {
            (line, Token::Keyword(k)) => Ok((line, k)),
            (line, t) => Err(Self::unexpected(line, "un mot-clé", &t)),
        }
    }

    fn keyword(&mut self, expected: &str, message: &str) -> Result<(), IoError> {
        match self.next()? {
            (_, Token::Keyword(k)) if k == expected => Ok(()),
            (line, _) => Err(error(format!("ligne {line} : {message}"))),
        }
    }

    /// `(valeur, valeur, …)`, la parenthèse ouvrante comprise.
    fn params(&mut self, line: usize) -> Result<Vec<Value>, IoError> {
        match self.next()? {
            (_, Token::Punct(b'(')) => {}
            (l, t) => return Err(Self::unexpected(l.max(line), "« ( »", &t)),
        }
        let mut out = Vec::new();
        if self.peek_punct(b')')? {
            self.next()?;
            return Ok(out);
        }
        loop {
            out.push(self.value()?);
            match self.next()? {
                (_, Token::Punct(b',')) => {}
                (_, Token::Punct(b')')) => return Ok(out),
                (l, t) => return Err(Self::unexpected(l, "« , » ou « ) »", &t)),
            }
        }
    }

    fn value(&mut self) -> Result<Value, IoError> {
        if self.peek_punct(b'(')? {
            let line = self.peek()?.0;
            return self.params(line).map(Value::List);
        }
        Ok(match self.next()? {
            (_, Token::Integer(i)) => Value::Integer(i),
            (_, Token::Real(x)) => Value::Real(x),
            (_, Token::String(s)) => Value::String(s),
            (_, Token::Enum(e)) => Value::Enum(e),
            (_, Token::Binary(b)) => Value::Binary(b),
            (_, Token::Instance(n)) => Value::Ref(n),
            (_, Token::Punct(b'$')) => Value::Null,
            (_, Token::Punct(b'*')) => Value::Derived,
            (line, Token::Keyword(k)) => Value::Typed(k, self.params(line)?),
            (line, t) => return Err(Self::unexpected(line, "une valeur", &t)),
        })
    }

    /// Instances jusqu'à `ENDSEC;`.
    fn data(&mut self, instances: &mut BTreeMap<u64, Instance>) -> Result<(), IoError> {
        loop {
            let (line, id) = match self.next()? {
                (_, Token::Keyword(k)) if k == "ENDSEC" => return self.punct(b';'),
                (line, Token::Instance(n)) => (line, n),
                (line, t) => {
                    return Err(Self::unexpected(
                        line,
                        "une instance « #n = » ou ENDSEC",
                        &t,
                    ))
                }
            };
            self.punct(b'=')?;
            let mut records = Vec::new();
            if self.peek_punct(b'(')? {
                // instance complexe : (NOM_A(…) NOM_B(…) …)
                self.next()?;
                while !self.peek_punct(b')')? {
                    let (l, name) = self.name()?;
                    records.push(Record {
                        params: self.params(l)?,
                        name,
                    });
                }
                self.next()?;
            } else {
                let (l, name) = self.name()?;
                records.push(Record {
                    params: self.params(l)?,
                    name,
                });
            }
            self.punct(b';')?;
            if instances.insert(id, Instance { line, records }).is_some() {
                return Err(error(format!(
                    "ligne {line} : l'instance #{id} est définie deux fois"
                )));
            }
        }
    }

    /// Saute une section facultative jusqu'à `ENDSEC;`.
    fn skip_section(&mut self) -> Result<(), IoError> {
        loop {
            match self.next()? {
                (_, Token::Keyword(k)) if k == "ENDSEC" => return self.punct(b';'),
                (line, Token::End) => {
                    return Err(error(format!("ligne {line} : fin du fichier avant ENDSEC")))
                }
                _ => {}
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cube unité, à la manière d'OCCT (en-tête AP214, instance complexe pour les unités, une face).
    const CUBE: &str = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('Open CASCADE Model'),'2;1');
FILE_NAME('cube','2026-10-28T10:00:00',('Tim'),(''),
  'Open CASCADE STEP processor 7.8','Open CASCADE 7.8','Unknown');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
/* une face plane du cube */
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = DIRECTION('',(0.,0.,1.));
#3 = DIRECTION('',(1.,0.,0.));
#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5 = PLANE('',#4);
#6 = ADVANCED_FACE('dessous',(#7),#5,.F.);
#7 = FACE_OUTER_BOUND('',#8,.T.);
#8 = EDGE_LOOP('',());
#9 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#10 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#9,
  'distance_accuracy_value','confusion accuracy');
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn cube_a_la_occt() {
        let p = parse(CUBE).unwrap();
        assert_eq!(p.schemas, ["AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }"]);
        assert_eq!(p.instances.len(), 10);
        let point = p.get(1, 0).unwrap();
        assert_eq!(point.line, 10);
        assert_eq!(
            point.records,
            [Record {
                name: "CARTESIAN_POINT".to_string(),
                params: vec![
                    Value::String(String::new()),
                    Value::List(vec![Value::Real(0.0); 3])
                ],
            }]
        );
        let face = &p.get(6, 0).unwrap().records[0];
        assert_eq!(face.params[0], Value::String("dessous".to_string()));
        assert_eq!(face.params[1], Value::List(vec![Value::Ref(7)]));
        assert_eq!(face.params[3], Value::Enum("F".to_string()));
        assert_eq!(
            p.get(8, 0).unwrap().records[0].params[1],
            Value::List(vec![])
        );
        // instance complexe : trois enregistrements, dans l'ordre du fichier
        let unit = p.get(9, 0).unwrap();
        let names: Vec<&str> = unit.records.iter().map(|r| r.name.as_str()).collect();
        assert_eq!(names, ["LENGTH_UNIT", "NAMED_UNIT", "SI_UNIT"]);
        assert_eq!(unit.record("NAMED_UNIT").unwrap().params, [Value::Derived]);
        assert_eq!(
            unit.record("SI_UNIT").unwrap().params,
            [Value::Enum("MILLI".into()), Value::Enum("METRE".into())]
        );
        // valeur typée sur deux lignes
        let uncertainty = p.get(10, 0).unwrap();
        assert_eq!(uncertainty.line, 19);
        assert_eq!(
            uncertainty.records[0].params[0],
            Value::Typed("LENGTH_MEASURE".into(), vec![Value::Real(1e-7)])
        );
    }

    fn data(body: &str) -> Result<Part21, IoError> {
        parse(&format!(
            "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n{body}\nENDSEC;\nEND-ISO-10303-21;\n"
        ))
    }

    fn params(body: &str) -> Vec<Value> {
        let p = data(&format!("#1 = X({body});")).unwrap();
        p.instances[&1].records[0].params.clone()
    }

    #[test]
    fn nombres() {
        assert_eq!(
            params("0, -12, +3, 1., -2.5E-03, 1.E+2, 6.02e23, 0.5"),
            [
                Value::Integer(0),
                Value::Integer(-12),
                Value::Integer(3),
                Value::Real(1.0),
                Value::Real(-2.5e-3),
                Value::Real(100.0),
                Value::Real(6.02e23),
                Value::Real(0.5),
            ]
        );
    }

    #[test]
    fn chaines_et_directives() {
        let s = |text: &str| match &params(text)[0] {
            Value::String(s) => s.clone(),
            v => panic!("{v:?}"),
        };
        assert_eq!(s("'l''axe'"), "l'axe");
        assert_eq!(s("'a\\\\b'"), "a\\b");
        assert_eq!(s("'\\X2\\00E9\\X0\\t\\X2\\00E900E8\\X0\\'"), "étéè");
        assert_eq!(s("'\\X\\E9'"), "é");
        assert_eq!(s("'\\S\\i'"), "é");
        assert_eq!(s("'\\X4\\0001F600\\X0\\'"), "😀");
        assert_eq!(s("'\\PA\\caf\\S\\i'"), "café");
        // une chaîne coupée sur deux lignes
        assert_eq!(s("'face\n_avant'"), "face_avant");
        // du texte UTF-8 brut passe tel quel
        assert_eq!(s("'pièce'"), "pièce");
    }

    #[test]
    fn autres_valeurs() {
        assert_eq!(
            params("$, *, .T., .UNSPECIFIED., #12, \"3F\", ((1, 2), ()), PARAMETER_VALUE(0.5)"),
            [
                Value::Null,
                Value::Derived,
                Value::Enum("T".into()),
                Value::Enum("UNSPECIFIED".into()),
                Value::Ref(12),
                Value::Binary("3F".into()),
                Value::List(vec![
                    Value::List(vec![Value::Integer(1), Value::Integer(2)]),
                    Value::List(vec![]),
                ]),
                Value::Typed("PARAMETER_VALUE".into(), vec![Value::Real(0.5)]),
            ]
        );
        assert_eq!(params(""), []);
    }

    #[test]
    fn plusieurs_sections_data_et_commentaires() {
        let text = "ISO-10303-21; /* en-tête */ HEADER; FILE_SCHEMA(('A','B')); ENDSEC;
DATA('partie 1',(('A'))); #1=X(1); ENDSEC;
DATA; #2 = /* au milieu */ Y(#1); ENDSEC;
END-ISO-10303-21;";
        let p = parse(text).unwrap();
        assert_eq!(p.schemas, ["A", "B"]);
        assert_eq!(p.instances.len(), 2);
        assert_eq!(p.instances[&2].records[0].params, [Value::Ref(1)]);
    }

    fn message(r: Result<Part21, IoError>) -> String {
        r.unwrap_err().to_string()
    }

    #[test]
    fn erreurs_avec_la_ligne() {
        let m = message(data("#1 = X(1)\n#2 = Y(2);"));
        assert!(m.contains("ligne 7") && m.contains("« ; » attendu"), "{m}");
        let m = message(data("#1 = X(1);\n#1 = Y(2);"));
        assert!(
            m.contains("ligne 7") && m.contains("#1 est définie deux fois"),
            "{m}"
        );
        let m = message(data("#1 = X('ouverte);"));
        assert!(
            m.contains("ligne 6") && m.contains("jamais refermée"),
            "{m}"
        );
        let m = message(data("#1 = X(1 2);"));
        assert!(m.contains("« , » ou « ) » attendu"), "{m}");
        let m = message(data("#1 = X(.T);"));
        assert!(m.contains("énumération"), "{m}");
        let m = message(data("#1 = X(1.2.3);"));
        assert!(m.contains("nombre mal formé « 1.2.3 »"), "{m}");
        let m = message(data("#1 = X(@);"));
        assert!(m.contains("caractère inattendu « @ »"), "{m}");
        let m = message(data("/* sans fin"));
        assert!(m.contains("commentaire"), "{m}");
        let m = message(data(""));
        assert!(m.contains("aucune instance"), "{m}");
        let p = data("#1 = X(#2);").unwrap();
        let m = p.get(2, 6).unwrap_err().to_string();
        assert!(m.contains("ligne 6") && m.contains("#2 est citée"), "{m}");
        assert!(m.starts_with("STEP : "), "{m}");
    }

    #[test]
    fn pas_un_step_texte() {
        assert!(message(parse("solid cube\nendsolid")).contains("ISO-10303-21"));
        assert!(message(parse("<?xml version=\"1.0\"?>")).contains("XML"));
        assert!(message(parse("PK\u{3}\u{4}")).contains("compressé"));
        let m = message(parse("ISO-10303-21;\nHEADER;\nENDSEC;\n"));
        assert!(m.contains("la fin du fichier"), "{m}");
    }
}
