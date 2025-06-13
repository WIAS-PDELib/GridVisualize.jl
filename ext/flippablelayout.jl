"""
This module is a submodule of `[GridVisualizeMakieExt](@ref)`.

It manages a layoutscene with interactive layout and blocking functionality.

Thanks to Julius Krumbiegel for providing  a basic implementation of focus switching.


`GridVisualize` avoids creating dependencies on plotting backends.
So we provide a way to emulate "import Makie" by allowing
to set it as a global variable in the [`GridVisualizeMakieExt.FlippableLayout.setmakie!`](@ref). 
As a consequence, we can't use Makie types at compile time.
"""
module FlippableLayout
using DocStringExtensions
using Observables

Makie = nothing

"""
$(TYPEDEF)


Struct describing flippable layout data.
We don't type annotate with Makie types as they are
unknown at start time.

$(FIELDS)
"""
mutable struct FLayout
    """
    Visible GridLayout
    """
    visible::Any #::GridLayout
    """
    Hidden GridLayot
    """
    offscreen::Any #::GridLayout

    """
    Has the layout been blocked by the block key ?
    """
    blocked::Bool

    """
    Layoutables attached to layout
    """
    layoutables::Dict{Tuple{Int64, Int64}, Any} # Union{Makie.Block,Makie.GridLayout}

    """
    Condition variable working together with the blocked field.
    """
    condition::Condition

    function FLayout(visible; blocked = false)
        return new(
            visible,
            Makie.GridLayout(; bbox = Makie.BBox(-500, -400, -500, -400)),
            blocked,
            Dict{Tuple{Int64, Int64}, Any}(),
            Condition()
        )
    end
end

function Base.setindex!(flayout::FLayout, layoutable, i, j)
    if isnothing(layoutable)
        flayout.offscreen[1, 1] = flayout.layoutables[(i, j)]
        delete!(flayout.layoutables, (i, j))
    elseif !isa(layoutable, Union{Makie.Block, Makie.GridLayout})
        error("can only set layoutables")
    else
        flayout.layoutables[(i, j)] = layoutable
        _showall(flayout)
        yield()
    end
    return flayout
end

Base.getindex(flayout::FLayout, i, j) = flayout.layoutables[(i, j)]

function _showall(flayout::FLayout)
    for (pos, layoutable) in flayout.layoutables
        flayout.visible[pos...] = layoutable
    end
    return Makie.trim!(flayout.visible)
end

#
# Check if mouse position is  within pixel area of scene
#
function _inarea(area, pos)
    return pos[1] > area.origin[1] &&
        pos[1] < area.origin[1] + area.widths[1] &&
        pos[2] > area.origin[2] &&
        pos[2] < area.origin[2] + area.widths[2]
end

function _inscene(l, pos)
    if isa(l, Makie.Block)
        return _inarea(l.scene.viewport[], pos)
    elseif isa(l, Makie.GridLayout)
        return _inarea(l.layoutobservables.computedbbox[], pos)
    end
end

"""
    flayoutscene(;blocked=false, kwargs...)

Layoutscene with interactive layout and blocking functionality.




The `,` key switches between focused view showing only one subscene
and "gallery view" showing all layoutables at once.

The space key toggles blocking of the execution of the main therad
when `yield` is replaced by `yieldwait`. Initial blocking state is 
set by the `blocked` kwarg.

The `kwargs...` are the same as of `AbstractPlotting.layoutscene`.

The idea is that this can work in some cases as a drop-in replacement
of `layoutscene`.     
"""
function flayoutscene(;
        blocked = false,
        focuskey = Makie.Keyboard.comma,
        blockingkey = Makie.Keyboard.space,
        kwargs...
    )
    parent = Makie.Scene(; camera = Makie.campixel!, kwargs...)
    layout = Makie.GridLayout(parent; alignmode = Makie.Outside(5))

    Makie.rowgap!(layout, Makie.Relative(0.0))
    Makie.colgap!(layout, Makie.Relative(0.0))

    flayout = FLayout(layout; blocked = blocked)

    gallery_view = Observable(true)

    # Watch mouse position
    mouseposition = Observable((0.0, 0.0))

    Makie.on(parent.events.mouseposition) do m
        mouseposition[] = m
        false
    end

    # Switch focus to subscene  at pos
    function _focus(focus)
        for (key, layoutable) in flayout.layoutables
            if key == focus
                flayout.visible[1, 1] = layoutable
            else
                flayout.offscreen[1, 1] = layoutable
            end
        end
        return Makie.trim!(flayout.visible)
    end

    # Figure out to which subscene the mouse position
    # corresponds
    function _subscene(mouseposition)
        for (key, layoutable) in flayout.layoutables
            if _inscene(layoutable, mouseposition)
                return key
            end
        end
        return nothing
    end

    function _toggle_block(flayout::FLayout)
        if flayout.blocked
            flayout.blocked = false
            notify(flayout.condition)
        else
            flayout.blocked = true
        end
        return nothing
    end

    # Handle global key events for `,` (focus/gallery view)
    # and space (toggle blocking)
    Makie.on(parent.events.keyboardbutton) do buttons
        if Makie.ispressed(parent, focuskey)
            if gallery_view[]
                pos = _subscene(mouseposition[])
                if !isnothing(pos)
                    _focus(pos)
                end
                gallery_view[] = false
            else
                _showall(flayout)
                gallery_view[] = true
            end
            return true
        end
        if Makie.ispressed(parent, blockingkey)
            _toggle_block(flayout)
            return true
        end
        return false
    end
    return parent, flayout
end

"""
     yieldwait(fliplayoutscene)

Yield and wait in case of scene being blocked via space key toggle
"""
function yieldwait(flayout::FLayout)
    yield()
    if flayout.blocked
        wait(flayout.condition)
    end
    return nothing
end

"""
    setmakie!(MyMakie)

Set the Makie module.
This Makie can be GLMakie,WGLMakie,CairoMakie
"""
setmakie!(MyMakie) = global Makie = MyMakie

end
