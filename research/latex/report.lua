-- Pandoc Lua filter for the PDF build:
-- 1. column widths of pipe tables proportional to their content
-- 2. TeX-like ^{...} and _{...} inside plain text become real super-/subscripts
local stringify = pandoc.utils.stringify

function Table(tbl)
  local ncols = #tbl.colspecs
  local maxlen = {}
  for i = 1, ncols do maxlen[i] = 0 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        local n = utf8.len(stringify(cell.contents)) or #stringify(cell.contents)
        if n > maxlen[i] then maxlen[i] = n end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do scan(body.body) end
  local total = 0
  for i = 1, ncols do total = total + maxlen[i] end
  if total < 70 then return tbl end          -- short table: natural widths
  local w, sum = {}, 0
  for i = 1, ncols do
    w[i] = math.max(math.min(maxlen[i], 55), 5) ^ 0.8
    sum = sum + w[i]
  end
  for i = 1, ncols do
    tbl.colspecs[i] = {tbl.colspecs[i][1], 0.97 * w[i] / sum}
  end
  return tbl
end

local function split(s)
  local out = {}
  local pos = 1
  while true do
    local a, b, kind, inner = s:find("([%^_]){([^{}]*)}", pos)
    if not a then break end
    if a > pos then table.insert(out, pandoc.Str(s:sub(pos, a - 1))) end
    if kind == "^" then
      table.insert(out, pandoc.Superscript(pandoc.Str(inner)))
    else
      table.insert(out, pandoc.Subscript(pandoc.Str(inner)))
    end
    pos = b + 1
  end
  if #out == 0 then return nil end
  if pos <= #s then table.insert(out, pandoc.Str(s:sub(pos))) end
  return out
end

-- x^y (single token) and single-letter bases with plain subscripts (σ_MC, ζ_G, x_v)
local function split_simple(s)
  -- superscript: base^token
  local a, b, pre, tok, post = s:find("^(.-%S)%^([%w%+%-−]+)(.*)$")
  if a and not s:find("://") then
    local out = {pandoc.Str(pre), pandoc.Superscript(pandoc.Str(tok))}
    if #post > 0 then table.insert(out, pandoc.Str(post)) end
    return out
  end
  -- subscript: exactly one (UTF-8) character before '_' at the start of the word
  local lead, base, tok2, post2 = s:match("^([%(|]?)([%z\1-\127\194-\244][\128-\191]*)_([%w%+%-]+)(.*)$")
  if base and utf8.len(base) == 1 and not s:find("://") then
    local out = {}
    if #lead > 0 then table.insert(out, pandoc.Str(lead)) end
    table.insert(out, pandoc.Str(base)); table.insert(out, pandoc.Subscript(pandoc.Str(tok2)))
    if #post2 > 0 then table.insert(out, pandoc.Str(post2)) end
    return out
  end
  return nil
end

function Str(el)
  if el.text:find("[%^_]{") then return split(el.text) end
  if el.text:find("[%^_]") then return split_simple(el.text) end
end
